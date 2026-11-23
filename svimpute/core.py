"""Compatible per-ALT ridge imputation for a jointly genotyped diploid locus.

This module does not repair incompatible reference genotypes. It consumes an
explicit, complete ALT catalogue, joint 0/1/2 copy counts, and the configured
outer/inner split specification. Missing calls remain missing in the truth.
No allele, sample, quality, frequency, or accuracy threshold is applied.

The configured seven ridge penalties are selected independently per ALT from
inner-fold SSE of predictions clipped to [0, 2]. Chosen ALT predictions are
then projected jointly onto {d >= 0, sum(d) <= 2}. Both forms are returned.
Projection supplies compatible expected copy counts, not genotype posterior
probabilities, hard calls, or evidence that the reference calls are correct.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from typing import Any

import numpy as np

ALPHAS = np.array([0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0])


def nearest_predictors(snp_chromosomes, snp_positions, chromosome, position,
                       width=64):
    """Same-chromosome nearest SNPs; stable distance ties, original order.

    All available SNPs are retained if the chromosome has fewer than width.
    This is the configured predictor-neighbourhood definition, not QC.
    """
    if width < 1:
        raise ValueError("Neighbour count must be positive")
    chromosomes = np.asarray(snp_chromosomes).astype(str)
    positions = np.asarray(snp_positions, dtype=np.int64)
    if chromosomes.shape != positions.shape or positions.ndim != 1:
        raise ValueError("SNP coordinate arrays must be aligned one-dimensional arrays")
    candidates = np.flatnonzero(chromosomes == str(chromosome))
    order = np.argsort(np.abs(positions[candidates] - int(position)), kind="stable")
    return np.sort(candidates[order[:width]])


def project_alt(dosages):
    """Euclidean projection on nonnegative ALT dosages with total <= 2.

    Any incomplete input row returns a wholly missing row. No hard genotype
    is inferred and more than two positive fractional dosages are permitted.
    """
    original = np.asarray(dosages, dtype=float)
    if original.ndim < 1 or original.shape[-1] == 0:
        raise ValueError("At least one ALT column is required")
    flat = original.reshape(-1, original.shape[-1])
    result = np.full_like(flat, np.nan)
    complete = np.isfinite(flat).all(axis=1)
    values = np.maximum(flat[complete], 0.0)
    excess = values.sum(axis=1) > 2.0
    for row in np.flatnonzero(excess):
        ordered = np.sort(values[row])[::-1]
        thresholds = (np.cumsum(ordered) - 2.0) / np.arange(1, len(ordered) + 1)
        rho = np.flatnonzero(ordered > thresholds)[-1]
        values[row] = np.maximum(values[row] - thresholds[rho], 0.0)
        # Remove at most floating-point excess, rather than exporting REF < 0.
        total = values[row].sum()
        if total > 2.0:
            values[row] *= 2.0 / total
    result[complete] = values
    return result.reshape(original.shape)


def reference_status(alt_counts):
    """Classify every reference sample without modifying its stored counts."""
    y = np.asarray(alt_counts, dtype=float)
    if y.ndim != 2 or not y.shape[1]:
        raise ValueError("ALT counts must have shape (reference samples, ALTs)")
    status = np.full(y.shape[0], "JOINT_GT_AVAILABLE", dtype=object)
    for i, row in enumerate(y):
        finite = np.isfinite(row)
        if np.isinf(row).any():
            status[i] = "INFINITE_REFERENCE_VALUE"
        elif not finite.any():
            status[i] = "MISSING_JOINT_GT"
        elif not finite.all():
            status[i] = "PARTIAL_JOINT_GT_UNCOMPUTABLE"
        elif (row < 0).any() or (row > 2).any():
            status[i] = "REFERENCE_COPY_COUNT_OUT_OF_RANGE"
        elif (row != np.floor(row)).any():
            status[i] = "REFERENCE_IS_NOT_HARD_JOINT_GT"
        elif row.sum() > 2:
            status[i] = "REFERENCE_EXCEEDS_DIPLOID_COPIES"
    return status


def validate_splits(split_specs, n_samples, n_outer=5, n_inner=3):
    """Require the configured nested explicit splits.

    An original sample-order check must precede passing positional indices.
    The public entry point additionally requires that order as an argument.
    """
    if n_outer < 2 or n_inner < 2 or len(split_specs) != n_outer * (n_inner + 1):
        raise ValueError("Expected n_outer x (n_inner + one outer) explicit splits")
    converted = []
    universe = np.arange(n_samples)
    for spec in split_specs:
        fit = np.asarray(spec["fit"], dtype=int)
        evaluate = np.asarray(spec["evaluate"], dtype=int)
        for values in (fit, evaluate):
            if values.ndim != 1 or len(np.unique(values)) != len(values):
                raise ValueError("Split indices must be unique one-dimensional arrays")
            if ((values < 0) | (values >= n_samples)).any():
                raise ValueError("Split index outside reference sample order")
        if np.intersect1d(fit, evaluate).size:
            raise ValueError("Training/evaluation overlap would leak labels")
        converted.append({"fit": fit, "evaluate": evaluate})
    outer_evaluation = []
    for outer in range(n_outer):
        final = converted[outer * (n_inner + 1) + n_inner]
        if not np.array_equal(np.sort(np.r_[final["fit"], final["evaluate"]]), universe):
            raise ValueError("Each outer split must partition all reference samples")
        inner_evaluation = []
        for inner in range(n_inner):
            spec = converted[outer * (n_inner + 1) + inner]
            if not np.array_equal(np.sort(np.r_[spec["fit"], spec["evaluate"]]),
                                  np.sort(final["fit"])):
                raise ValueError("Inner split must partition only its outer training set")
            inner_evaluation.extend(spec["evaluate"].tolist())
        if not np.array_equal(np.sort(inner_evaluation), np.sort(final["fit"])):
            raise ValueError("Inner evaluation must visit each outer-training sample once")
        outer_evaluation.extend(final["evaluate"].tolist())
    if not np.array_equal(np.sort(outer_evaluation), universe):
        raise ValueError("Outer evaluation must visit each reference sample once")
    return converted


def _scale(train_x, evaluate_x):
    finite = np.isfinite(train_x)
    counts = finite.sum(axis=0)
    mean = np.divide(np.where(finite, train_x, 0.0).sum(axis=0), counts,
                     out=np.zeros(train_x.shape[1]), where=counts > 0)
    filled = np.where(finite, train_x, mean)
    sd = filled.std(axis=0)
    constant = sd == 0
    sd[constant] = 1.0
    zfit = (filled - mean) / sd
    zevaluate = (np.where(np.isfinite(evaluate_x), evaluate_x, mean) - mean) / sd
    return zfit, zevaluate, {"mean": mean, "sd": sd,
                            "observed_predictor_count": counts,
                            "constant_predictor": constant}


def _curve(x, y, fit, evaluate_x, available, alphas):
    """Predict all penalties with scaling fitted on the split's SNP rows only."""
    result = np.full((len(alphas), len(evaluate_x), y.shape[1]), np.nan)
    zfit, zevaluate, scaling = _scale(x[fit], evaluate_x)
    fit_mask = available[fit]
    diagnostics = {"n_training_samples": int(len(fit)),
                   "n_joint_gt_training": int(fit_mask.sum()),
                   "all_missing_predictors": np.flatnonzero(
                       scaling["observed_predictor_count"] == 0).tolist()}
    if not fit_mask.any():
        diagnostics["status"] = "NO_JOINT_GT_IN_TRAINING"
        return result, diagnostics, scaling, None
    xx, yy = zfit[fit_mask], y[fit[fit_mask]]
    xmean, ymean = xx.mean(axis=0), yy.mean(axis=0)
    xc, yc = xx - xmean, yy - ymean
    # Thin SVD avoids a p-by-p matrix when the user supplies a large SNP pool.
    _, singular, vt = np.linalg.svd(xc, full_matrices=False)
    eigenvalues, eigenvectors = singular ** 2, vt.T
    left = (zevaluate - xmean) @ eigenvectors
    right = eigenvectors.T @ (xc.T @ yc)
    inv = 1.0 / (eigenvalues[:, None] + alphas[None, :])
    result = np.einsum("ei,ia,im->aem", left, inv, right, optimize=True) + ymean
    diagnostics["status"] = "COMPUTED"
    parameters = {"eigenvectors": eigenvectors, "right": right,
                  "eigenvalues": eigenvalues, "xmean": xmean, "ymean": ymean}
    return np.clip(result, 0.0, 2.0), diagnostics, scaling, parameters


def allele_metrics(truth, prediction, allele_ids):
    """Signed Pearson dosage correlation and MSE; no selection threshold."""
    truth, prediction = np.asarray(truth), np.asarray(prediction)
    records = []
    for j, allele in enumerate(allele_ids):
        ok = np.isfinite(truth[:, j]) & np.isfinite(prediction[:, j])
        a, b = truth[ok, j], prediction[ok, j]
        record = {"allele_id": str(allele), "n_pairs": int(ok.sum()),
                  "n_truth_available": int(np.isfinite(truth[:, j]).sum()),
                  "pearson_r": None, "mse": None}
        if not len(a):
            record["pearson_status"] = "NO_PAIRED_VALUES"
        else:
            record["mse"] = float(np.mean((a - b) ** 2))
            if len(a) < 2:
                record["pearson_status"] = "FEWER_THAN_TWO_PAIRED_VALUES"
            elif np.ptp(a) == 0:
                record["pearson_status"] = "CONSTANT_JOINT_GT"
            elif np.ptp(b) == 0:
                record["pearson_status"] = "CONSTANT_PREDICTED_DOSAGE"
            else:
                record["pearson_r"] = float(np.corrcoef(a, b)[0, 1])
                record["pearson_status"] = "COMPUTED"
        records.append(record)
    return records


def compatibility_metrics(dosages):
    dosages = np.asarray(dosages)
    complete = np.isfinite(dosages).all(axis=1)
    values = dosages[complete]
    return {"n_rows": len(dosages), "n_complete": int(complete.sum()),
            "n_uncomputable": int((~complete).sum()),
            "n_negative_alt": int((values < 0).any(axis=1).sum()),
            "n_sum_alt_exceeds_two": int((values.sum(axis=1) > 2.0).sum()),
            "max_sum_alt": float(values.sum(axis=1).max()) if len(values) else None}


def impute_locus(*, x_reference, x_target, alt_counts, split_specs,
                 reference_ids, split_sample_ids, target_ids, allele_ids,
                 feature_indices, snp_ids, locus_id, reference_provenance,
                 alphas=ALPHAS, n_outer=5, n_inner=3):
    """Return nested OOF and target predictions for one complete ALT catalogue.

    `reference_provenance` must describe the new joint genotyping source. The
    caller must not pass the old independent per-record GT as a corrected truth.
    IDs, original truth, selected predictors, all preprocessing, per-ALT alpha,
    fold diagnostics, and fitted coefficients are retained in the result.
    """
    x = np.asarray(x_reference, dtype=float)
    xt = np.asarray(x_target, dtype=float)
    y = np.asarray(alt_counts, dtype=float)
    indices = np.asarray(feature_indices, dtype=int)
    alphas = np.asarray(alphas, dtype=float)
    reference_ids = [str(i) for i in reference_ids]
    split_sample_ids = [str(i) for i in split_sample_ids]
    target_ids = [str(i) for i in target_ids]
    allele_ids = [str(i) for i in allele_ids]
    if reference_ids != split_sample_ids:
        raise ValueError("Reference sample order differs from configured split order")
    for ids, name in [(reference_ids, "reference"), (target_ids, "target"),
                      (allele_ids, "allele")]:
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate {name} IDs require explicit resolution")
    if not reference_provenance:
        raise ValueError("Joint reference genotype provenance is required")
    if x.ndim != 2 or xt.ndim != 2 or x.shape[1] != xt.shape[1]:
        raise ValueError("SNP matrices must share aligned predictor columns")
    if x.shape[0] != len(reference_ids) or xt.shape[0] != len(target_ids):
        raise ValueError("SNP sample IDs and matrix rows do not align")
    if y.shape != (len(reference_ids), len(allele_ids)):
        raise ValueError("ALT count matrix does not align with sample/allele IDs")
    if len(snp_ids) != x.shape[1]:
        raise ValueError("SNP IDs and matrix columns do not align")
    if indices.ndim != 1 or len(indices) != len(np.unique(indices)):
        raise ValueError("Predictor indices must be unique")
    if ((indices < 0) | (indices >= x.shape[1])).any():
        raise ValueError("Predictor index outside SNP matrix")
    if alphas.ndim != 1 or not len(alphas) or not np.isfinite(alphas).all() or (alphas <= 0).any():
        raise ValueError("Ridge penalties must be finite and positive")
    specs = validate_splits(split_specs, len(reference_ids), n_outer, n_inner)
    x, xt = x[:, indices], xt[:, indices]
    rows = reference_status(y)
    available = rows == "JOINT_GT_AVAILABLE"
    incomplete = np.isin(rows, ["MISSING_JOINT_GT", "PARTIAL_JOINT_GT_UNCOMPUTABLE"])
    result: dict[str, Any] = {
        "locus_id": str(locus_id), "reference_provenance": reference_provenance,
        "reference_ids": reference_ids, "target_ids": target_ids,
        "allele_ids": allele_ids, "truth": y.copy(), "reference_row_status": rows,
        "reference_status_counts": dict(Counter(rows.tolist())),
        "feature_indices": indices, "feature_snp_ids": np.asarray(snp_ids)[indices],
        "n_predictors": len(indices), "alphas": alphas.copy(),
        "oof_unconstrained": np.full_like(y, np.nan),
        "oof_compatible": np.full_like(y, np.nan),
        "target_unconstrained": np.full((len(xt), y.shape[1]), np.nan),
        "target_compatible": np.full((len(xt), y.shape[1]), np.nan),
        "outer_alpha_indices": np.full((n_outer, y.shape[1]), -1, dtype=int),
        "final_alpha_indices": np.full(y.shape[1], -1, dtype=int),
        "inner_sse": np.zeros((n_outer, len(alphas), y.shape[1])),
        "inner_scored_counts": np.zeros(n_outer, dtype=int),
        "split_diagnostics": [], "preprocessing": [],
        "no_accuracy_selection": True, "no_quality_exclusions": True,
        "alpha_selection": "per ALT pooled inner SSE of clipped predictions",
        "constraint": "Euclidean projection after configured clipping; nonnegative ALT sum <= 2",
    }
    if not (available | incomplete).all():
        result["status"] = "REFERENCE_INCOMPATIBLE_JOINT_GT_NOT_FITTED"
    elif len(indices) == 0:
        result["status"] = "NO_SAME_CHROMOSOME_SNP_PREDICTORS"
    elif not available.any():
        result["status"] = "NO_COMPLETE_JOINT_REFERENCE_GT"
    else:
        for outer in range(n_outer):
            ss = result["inner_sse"][outer]
            for inner in range(n_inner):
                number = outer * (n_inner + 1) + inner
                spec = specs[number]
                ev = spec["evaluate"]
                pred, diag, scaling, _ = _curve(
                    x, y, spec["fit"], x[ev], available, alphas)
                score = available[ev] & np.isfinite(pred).all(axis=(0, 2))
                count = int(score.sum())
                if count:
                    ss += ((pred[:, score] - y[ev[score]]) ** 2).sum(axis=1)
                result["inner_scored_counts"][outer] += count
                diag.update(split_index=number, outer=outer, inner=inner,
                            n_joint_gt_evaluation=int(available[ev].sum()),
                            n_scored=count, fit=spec["fit"], evaluate=ev)
                result["split_diagnostics"].append(diag)
                result["preprocessing"].append(scaling)
            number = outer * (n_inner + 1) + n_inner
            spec = specs[number]
            ev = spec["evaluate"]
            pred, diag, scaling, _ = _curve(x, y, spec["fit"], x[ev], available, alphas)
            if result["inner_scored_counts"][outer] > 0 and np.isfinite(pred).all():
                alpha_index = ss.argmin(axis=0)
                chosen = pred[alpha_index, :, np.arange(y.shape[1])].T
                result["outer_alpha_indices"][outer] = alpha_index
                result["oof_unconstrained"][ev] = chosen
                result["oof_compatible"][ev] = project_alt(chosen)
            elif result["inner_scored_counts"][outer] == 0:
                diag["status"] = "NO_COMPUTABLE_INNER_VALIDATION"
            diag.update(split_index=number, outer=outer, inner=None,
                        fit=spec["fit"], evaluate=ev)
            result["split_diagnostics"].append(diag)
            result["preprocessing"].append(scaling)
        if result["inner_scored_counts"].sum() > 0:
            final_alpha = result["inner_sse"].sum(axis=0).argmin(axis=0)
            pred, diag, scaling, parameters = _curve(
                x, y, np.arange(len(x)), xt, available, alphas)
            chosen = pred[final_alpha, :, np.arange(y.shape[1])].T
            result["target_unconstrained"] = chosen
            result["target_compatible"] = project_alt(chosen)
            result["final_alpha_indices"] = final_alpha
            result["final_preprocessing"] = scaling
            result["final_diagnostics"] = diag
            if parameters is not None:
                denominator = parameters["eigenvalues"][:, None] + alphas[final_alpha][None, :]
                coefficient = parameters["eigenvectors"] @ (parameters["right"] / denominator)
                result["final_coefficient"] = coefficient
                result["final_intercept"] = parameters["ymean"] - parameters["xmean"] @ coefficient
            result["status"] = "COMPUTED" if np.isfinite(result["oof_compatible"]).all() else "PARTIALLY_COMPUTABLE_OOF"
        else:
            result["status"] = "NO_COMPUTABLE_INNER_VALIDATION"
    # Partial calls have no interpretable joint reference vector for scoring.
    scoring_truth = np.where(available[:, None], y, np.nan)
    for scope in ("oof", "target"):
        for kind in ("unconstrained", "compatible"):
            key = f"{scope}_{kind}"
            result[key + "_compatibility"] = compatibility_metrics(result[key])
            if scope == "oof":
                result[key + "_metrics"] = allele_metrics(scoring_truth, result[key], allele_ids)
        result[scope + "_ref_dosage"] = 2.0 - result[scope + "_compatible"].sum(axis=1)
    return result


def write_result(result, prefix):
    """Store all arrays losslessly and JSON metadata, with source hashes."""
    from pathlib import Path
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    arrays = {}

    def encode(value, key):
        if isinstance(value, np.ndarray):
            if value.dtype == object:
                value = value.astype(str)
            arrays[key] = value
            return {"array": key, "shape": list(value.shape), "dtype": str(value.dtype)}
        if isinstance(value, dict):
            return {k: encode(v, key + "__" + k) for k, v in value.items()}
        if isinstance(value, (tuple, list)):
            return [encode(v, key + "__" + str(i)) for i, v in enumerate(value)]
        if isinstance(value, np.generic):
            return value.item()
        return value

    metadata = encode(result, "result")
    array_path = prefix.with_suffix(".npz")
    np.savez_compressed(array_path, **arrays)
    metadata["array_file"] = array_path.name
    metadata["array_sha256"] = hashlib.sha256(array_path.read_bytes()).hexdigest()
    metadata["module_sha256"] = hashlib.sha256(__import__("pathlib").Path(__file__).read_bytes()).hexdigest()
    prefix.with_suffix(".json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return metadata
