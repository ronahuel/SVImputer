"""Choose a local SNP neighborhood from a mixed reference VCF and an LD panel."""
from collections import Counter, deque
from contextlib import ExitStack
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from .cli import make_splits, write_tsv
from .core import ALPHAS, allele_metrics, impute_locus, validate_splits
from .panel_vcf import match_panel, neighbors, target_loci
from .profile import ReferenceProfile
from .report import Accuracy, plot, writer


def folds(args, sample_ids):
    if args.splits_json:
        data = json.loads(Path(args.splits_json).read_text())
        if data["sample_ids"] != sample_ids or len(data["repeats"]) != args.repeats:
            raise ValueError("Split sample order or repeat count differs from the reference/options.")
        repeats = data["repeats"]
    else:
        repeats = [make_splits(len(sample_ids), args.outer_folds, args.inner_folds, args.seed + i)
                   for i in range(args.repeats)]
    for splits in repeats:
        validate_splits(splits, len(sample_ids), args.outer_folds, args.inner_folds)
    return repeats


def dosage_row(values):
    return [format(v, ".8g") if np.isfinite(v) else "." for v in values]


def evaluate(args, out):
    panel_format = "bim" if args.ld_panel_bim else "vcf"
    panel_path = args.ld_panel_bim or args.ld_panel_vcf
    profile = ReferenceProfile(out / "reference_profile")
    names, new_names, markers, usable, counts = match_panel(args.reference, panel_path, args.chromosome,
                                                            panel_format, profile=profile)
    profile.save()
    panels = Counter(m.chrom for m in markers)
    chroms = sorted(set(counts) | set(panels) | set(args.chromosome or []))
    split_sets = folds(args, names)
    (out / "splits.json").write_text(json.dumps({"sample_ids": names, "repeats": split_sets}))
    if panel_format == "bim":
        write_tsv(out / "ld_matches.tsv",
                  (dict(chrom=m.chrom, pos=m.pos, id=m.marker_id, a1=m.alt, a2=m.ref, status=m.status,
                        reference_ref=m.reference_pair[0] if m.status == "MATCHED" else "",
                        reference_alt=m.reference_pair[1] if m.status == "MATCHED" else "") for m in markers),
                  ["chrom", "pos", "id", "a1", "a2", "status", "reference_ref", "reference_alt"])
    else:
        write_tsv(out / "ld_matches.tsv",
                  (dict(chrom=m.chrom, pos=m.pos, ref=m.ref, alt=m.alt, status=m.status,
                        reversed_alleles=m.flipped) for m in markers),
                  ["chrom", "pos", "ref", "alt", "status", "reversed_alleles"])
    maximum = max(map(len, usable.values()), default=0)
    # Above this maximum all chromosomes use exactly the same features.
    grid = sorted({min(k, maximum) if maximum else k for k in args.neighbors})
    coords = {chrom: np.array([m.pos for m in group]) for chrom, group in usable.items()}
    panel_sites = {(m.chrom, m.pos, m.reference_pair[0], (m.reference_pair[1],))
                   for group in usable.values() for m in group}
    target_counts = Counter()
    acc = Accuracy()
    with ExitStack() as stack:
        def table(filename, fields, compressed=False):
            handle = (gzip.open(out / filename, "wt", newline="") if compressed
                      else (out / filename).open("w", newline=""))
            return writer(stack.enter_context(handle), fields)

        metrics = table("accuracy.tsv.gz", ["chrom", "pos", "id", "ref", "alt", "repeat", "k",
                         "used", "pairs", "pearson_r", "reason"], True)
        feature_table = table("neighbors.tsv.gz", ["chrom", "pos", "ld_positions_nearest_first"], True)
        prefix = ["chrom", "pos", "ref", "alt", "repeat", "k"]
        oof = table("dosage_cv.tsv.gz", prefix + names, True)
        new = table("dosage_ld.tsv.gz", prefix + new_names, True) if new_names else None
        for locus in target_loci(args.reference, names, panel_sites, args.chromosome, args.targets):
            chrom, pos = locus["chrom"], locus["pos"]
            target_counts[chrom] += 1
            ranked = neighbors(usable.get(chrom, []), coords.get(chrom, []), pos, max(grid))
            feature_table.writerow([chrom, pos, ",".join(str(m.pos) for m in ranked)])
            sizes = {k: min(k, len(ranked)) for k in grid}
            correlations = {k: {} for k in grid}
            alt_ids = locus["alts"]
            source_ids = ",".join(r["id"] for r in locus["source"])
            ref = locus["ref"]
            for used in sorted(set(sizes.values())):
                caps = [k for k in grid if sizes[k] == used]
                reason = locus["reason"] or (None if used else "NO_LOCAL_LD_SNPS")
                if not reason:
                    features = ranked[:used]
                    x = np.column_stack([m.ref_calls for m in features])
                    xt = np.column_stack([m.new_calls for m in features])
                for repeat, splits in enumerate(split_sets, 1):
                    if reason:
                        for k in caps:
                            for alt in alt_ids or ["."]:
                                metrics.writerow([chrom, pos, source_ids, ref, alt, repeat, k, used, 0, "", reason])
                        continue
                    res = impute_locus(x_reference=x, x_target=xt, alt_counts=locus["counts"],
                        split_specs=splits, reference_ids=names, split_sample_ids=names, target_ids=new_names,
                        allele_ids=alt_ids, feature_indices=np.arange(used), snp_ids=[str(m.pos) for m in features],
                        locus_id=locus["locus"], reference_provenance={"GT": "mixed reference VCF"},
                        alphas=ALPHAS, n_outer=args.outer_folds, n_inner=args.inner_folds)
                    # A sample without any local SNP call has no basis for an individual prediction.
                    held, predicted = res["oof_compatible"], res["target_compatible"]
                    held[~np.isfinite(x).any(axis=1)] = np.nan
                    predicted[~np.isfinite(xt).any(axis=1)] = np.nan
                    complete = np.isfinite(locus["counts"]).all(axis=1)
                    truth = np.where(complete[:, None], locus["counts"], np.nan)
                    scores = allele_metrics(truth, held, alt_ids)
                    for k in caps:
                        for j, score in enumerate(scores):
                            value = score["pearson_r"]
                            state = score["pearson_status"]
                            if state == "NO_PAIRED_VALUES" and res["status"] != "COMPUTED":
                                state = res["status"]
                            alt = alt_ids[j]
                            metrics.writerow([chrom, pos, source_ids, ref, alt, repeat, k, used,
                                              score["n_pairs"], value, state])
                            if value is not None:
                                correlations[k][repeat, j] = value
                            key = [chrom, pos, ref, alt, repeat, k]
                            oof.writerow(key + dosage_row(held[:, j]))
                            if new:
                                new.writerow(key + dosage_row(predicted[:, j]))
            acc.add_locus(chrom, correlations, sizes)
    fields = ["chrom", "reference_snps", "reference_sv", "reference_indels", "other_records",
              "ld_supplied", "ld_matched", "target_loci"]
    write_tsv(out / "chromosomes.tsv",
              (dict(zip(fields, [c, counts[c]["snp"], counts[c]["sv"], counts[c]["indel"],
                                 counts[c]["other"] + counts[c]["unresolved"], panels[c],
                                 len(usable.get(c, [])), target_counts[c]])) for c in chroms), fields)
    rows = acc.rows()
    fields = ["scope", "chrom", "k", "min_used", "max_used", "loci", "evaluable_loci",
              "paired_loci", "mean_r_available", "mean_r_paired", "best"]
    write_tsv(out / "window_accuracy.tsv", rows, fields)
    best = [r["k"] for r in rows if r["scope"] == "genome" and r["best"]]
    plot(rows, out / "accuracy_by_neighbors.png")
    report = dict(reference=args.reference, ld=panel_path, ld_panel_format=panel_format,
                  reference_samples=len(names), new_samples=len(new_names),
                  targets=args.targets, requested_neighbors=args.neighbors, evaluated_caps=grid,
                  outer_folds=args.outer_folds, inner_folds=args.inner_folds, repeats=args.repeats,
                  seed=args.seed, alphas=ALPHAS.tolist(), chromosomes=args.chromosome, threads=args.threads,
                  target_loci=sum(target_counts.values()), matched_ld_snps=sum(map(len, usable.values())),
                  best_caps=best, state="COMPLETE" if best else "NO_COMPARABLE_ACCURACY")
    (out / "run.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def run(args, *, describe_only=False):
    out = Path(args.out)
    if out.exists():
        raise ValueError(f"Choose a new output directory. Already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=out.parent, prefix=".sv-impute-") as staging:
        payload = Path(staging) / "results"
        payload.mkdir()
        if describe_only:
            import pysam
            profile = ReferenceProfile(payload)
            with pysam.VariantFile(args.reference) as source:
                deque(profile.scan(source, args.chromosome), maxlen=0)
            report = profile.save()
        else:
            report = evaluate(args, payload)
        payload.rename(out)
    if describe_only:
        print(f"Reference profile: {out}")
        return report
    print(f"{report['target_loci']} target loci; {report['matched_ld_snps']} matched LD SNPs.")
    best = ", ".join(map(str, report["best_caps"])) or "undefined"
    print(f"Best neighborhood cap(s): {best}. Results: {out}")
    return report
