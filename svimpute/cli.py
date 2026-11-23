"""Species-independent chromosome labels; diploid GT input only, no QC filtering."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path


def checksum(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_list(path):
    if path is None:
        return None
    rows = [line.strip() for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    if len(rows) != len(set(rows)):
        raise ValueError('Duplicate entries in list: ' + str(path))
    return rows


def load_vcf(path, pool=None, chromosomes=None, targets=None):
    import numpy as np
    import pysam
    records, matched, seen_keys = [], Counter(), set()
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        if not samples or len(samples) != len(set(samples)):
            raise ValueError('VCF needs unique sample IDs: ' + str(path))
        for rec in vcf:
            if chromosomes and rec.contig not in chromosomes:
                continue
            alts = list(rec.alts or [])
            key = f'{rec.contig}:{rec.pos}:{rec.ref}:{",".join(alts)}'
            selectors = {key}
            if rec.id and rec.id != '.':
                selectors.add(rec.id)
            requested = pool if pool is not None else targets
            if requested is not None and not selectors.intersection(requested):
                continue
            for name in selectors.intersection(requested or []):
                matched[name] += 1
            if key in seen_keys:
                raise ValueError('Duplicate VCF record requires resolution: ' + key)
            seen_keys.add(key)
            counts = np.full((len(samples), len(alts)), np.nan)
            statuses, calls = [], []
            for i, sample in enumerate(samples):
                gt = rec.samples[sample].get('GT')
                calls.append(list(gt) if gt is not None else None)
                if gt is None or all(a is None for a in gt):
                    status = 'MISSING_GT'
                elif len(gt) != 2:
                    status = 'NON_DIPLOID_GT_UNSUPPORTED'
                elif any(a is None for a in gt):
                    status = 'PARTIAL_GT_UNCOMPUTABLE'
                elif any(a < 0 or a > len(alts) for a in gt):
                    status = 'GT_ALLELE_INDEX_OUT_OF_RANGE'
                else:
                    status = 'CALLED'
                    counts[i] = [sum(a == j for a in gt) for j in range(1, len(alts) + 1)]
                statuses.append(status)
            records.append({'key': key, 'id': rec.id, 'chrom': rec.contig,
                            'pos': rec.pos, 'ref': rec.ref, 'alts': alts,
                            'end': rec.stop, 'counts': counts,
                            'sample_status': statuses, 'original_GT': calls})
    if pool is not None or targets is not None:
        requested = pool if pool is not None else targets
        failures = {name: matched[name] for name in requested if matched[name] != 1}
        if failures:
            raise ValueError('Requested IDs absent or ambiguous in this chromosome scope: ' + str(failures))
    return samples, records


def make_splits(n, outer, inner, seed):
    import numpy as np
    if outer < 2 or inner < 2 or n < outer or n - (n + outer - 1) // outer < inner:
        raise ValueError('Sample count cannot support the requested nested folds')
    rng = np.random.default_rng(seed)
    universe = np.arange(n)
    result = []
    for evaluation in np.array_split(rng.permutation(n), outer):
        fit = np.setdiff1d(universe, evaluation)
        for inner_eval in np.array_split(rng.permutation(fit), inner):
            result.append({'fit': np.setdiff1d(fit, inner_eval).tolist(), 'evaluate': inner_eval.tolist()})
        result.append({'fit': fit.tolist(), 'evaluate': evaluation.tolist()})
    return result


def locus_records(records):
    import numpy as np
    groups = defaultdict(list)
    for record in records:
        groups[(record['chrom'], record['pos'])].append(record)
    for (chrom, pos), rows in groups.items():
        alleles = [a for r in rows for a in r['alts']]
        reason = None
        if not alleles:
            reason = 'NO_ALT_ALLELE'
        elif len(set(alleles)) != len(alleles):
            reason = 'DUPLICATED_ALT_REPRESENTATION_UNRESOLVED'
        elif len({(r['ref'], r['end']) for r in rows}) != 1:
            reason = 'POSITIONAL_RECORDS_REQUIRE_COMMON_REFERENCE_REPRESENTATION'
        elif any(s in {'NON_DIPLOID_GT_UNSUPPORTED', 'GT_ALLELE_INDEX_OUT_OF_RANGE'}
                 for r in rows for s in r['sample_status']):
            reason = 'UNSUPPORTED_REFERENCE_GT_FORMAT'
        counts = np.concatenate([r['counts'] for r in rows], axis=1)
        if np.any(np.nansum(counts, axis=1) > 2):
            reason = 'REFERENCE_EXCEEDS_DIPLOID_COPIES'
        yield {'locus': f'{chrom}:{pos}', 'chrom': chrom, 'pos': pos,
               'ref': rows[0]['ref'], 'end': max(r['end'] for r in rows),
               'alts': alleles, 'source': rows, 'reason': reason,
               'counts': counts}


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    for key in ['OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS']:
        os.environ[key] = str(args.threads)
    import numpy as np
    from .core import ALPHAS, impute_locus, nearest_predictors, write_result
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Output directory must be empty; source files are never overwritten')
    out.mkdir(parents=True, exist_ok=True)
    pool, wanted = read_list(args.predictor_pool), read_list(args.target_ids)
    # Pool is loaded across chromosomes first, so a genome-wide panel remains valid
    # when only one target chromosome is requested.
    snp_samples, predictors = load_vcf(args.snps_vcf, pool=pool)
    ref_samples, targets = load_vcf(args.targets_vcf, chromosomes=args.chromosome, targets=wanted)
    if not targets:
        raise ValueError('No target variants in the explicitly requested scope')
    snp_lookup = {sample: i for i, sample in enumerate(snp_samples)}
    missing = [s for s in ref_samples if s not in snp_lookup]
    if missing:
        raise ValueError('Target-reference samples absent from predictor VCF: ' + ','.join(missing))
    predictor_errors = []
    for r in predictors:
        if len(r['ref']) != 1 or len(r['alts']) != 1 or len(r['alts'][0]) != 1 or r['ref'] not in 'ACGT' or r['alts'][0] not in 'ACGT':
            predictor_errors.append({'variant': r['key'], 'reason': 'PREDICTOR_REQUIRES_BIALLELIC_SNP'})
        if any(s in {'NON_DIPLOID_GT_UNSUPPORTED', 'GT_ALLELE_INDEX_OUT_OF_RANGE'} for s in r['sample_status']):
            predictor_errors.append({'variant': r['key'], 'reason': 'UNSUPPORTED_PREDICTOR_GT'})
    if predictor_errors:
        (out/'input_errors.json').write_text(json.dumps(predictor_errors, indent=2))
        raise ValueError('Predictor formats require resolution; see input_errors.json. No records silently removed.')
    target_samples = [s for s in snp_samples if s not in set(ref_samples)]
    matrix = np.column_stack([r['counts'][:, 0] for r in predictors]) if predictors else np.empty((len(snp_samples), 0))
    x = matrix[[snp_lookup[s] for s in ref_samples]]
    xt = matrix[[snp_lookup[s] for s in target_samples]]
    chromosomes = np.array([r['chrom'] for r in predictors])
    positions = np.array([r['pos'] for r in predictors], dtype=np.int64)
    predictor_ids = [r['key'] for r in predictors]
    if args.splits_json:
        split_doc = json.loads(Path(args.splits_json).read_text())
        if split_doc['sample_ids'] != ref_samples:
            raise ValueError('Supplied split sample order must match target VCF order exactly')
        repeats = split_doc['repeats']
        if len(repeats) != args.repeats:
            raise ValueError('Supplied split repeat count differs from --repeats')
    else:
        repeats = [make_splits(len(ref_samples), args.outer_folds, args.inner_folds, args.seed + rep) for rep in range(args.repeats)]
    (out/'splits.json').write_text(json.dumps({'sample_ids': ref_samples, 'repeats': repeats}, indent=2))
    metrics, statuses, chosen_features = [], [], []
    ds_paths = []
    for rep in range(args.repeats):
        ds_paths.append(out/f'imputed_repeat{rep+1:02d}.vcf')
    streams = [path.open('w') for path in ds_paths] if not args.check_only and target_samples else []
    for stream in streams:
        stream.write('##fileformat=VCFv4.2\n##source=sv-impute-framework-0.1.0\n')
        for chrom in dict.fromkeys(r['chrom'] for r in targets):
            stream.write(f'##contig=<ID={chrom}>\n')
        stream.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Source reference end">\n')
        stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Not inferred from expected dosage">\n')
        stream.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="Compatible expected alternate allele copy counts">\n')
        stream.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' + ''.join('\t'+s for s in target_samples) + '\n')
    try:
        for number, locus in enumerate(locus_records(targets), 1):
            if args.predictor_mode == 'nearest':
                features = nearest_predictors(chromosomes, positions, locus['chrom'], locus['pos'], args.neighbors)
            elif args.predictor_mode == 'chromosome':
                features = np.flatnonzero(chromosomes == locus['chrom'])
            else:
                features = np.arange(len(predictors))
            reason = locus['reason']
            # Refuse self-prediction instead of silently changing the user's pool.
            if any(predictors[j]['chrom'] == locus['chrom'] and predictors[j]['pos'] == locus['pos'] for j in features):
                reason = 'TARGET_POSITION_IN_PREDICTORS_RESOLVE_POOL_TO_PREVENT_LEAKAGE'
            chosen_features.append({'locus': locus['locus'], 'n_predictors': len(features), 'predictors': ';'.join(predictor_ids[j] for j in features)})
            input_doc = {'locus': locus['locus'], 'sample_ids': ref_samples,
                         'source_records': [{k: v for k, v in r.items() if k != 'counts'} for r in locus['source']],
                         'input_reason': reason}
            folder = out/'loci'/f'{number:08d}';folder.mkdir(parents=True)
            (folder/'input.json').write_text(json.dumps(input_doc, indent=2))
            for rep, splits in enumerate(repeats):
                result = None
                state = reason or 'INPUT_READY'
                if not args.check_only and reason is None:
                    result = impute_locus(x_reference=x, x_target=xt, alt_counts=locus['counts'],
                        split_specs=splits, reference_ids=ref_samples, split_sample_ids=ref_samples,
                        target_ids=target_samples, allele_ids=[f'ALT{i+1}' for i in range(len(locus['alts']))],
                        feature_indices=features, snp_ids=predictor_ids, locus_id=locus['locus'],
                        reference_provenance={'called_vcf': str(args.targets_vcf), 'format': 'GT', 'positional_grouping': True},
                        alphas=np.array(args.alphas) if args.alphas else ALPHAS,
                        n_outer=args.outer_folds, n_inner=args.inner_folds)
                    write_result(result, folder/f'repeat{rep+1:02d}')
                    state = result['status']
                    for kind in ['unconstrained', 'compatible']:
                        for metric in result['oof_'+kind+'_metrics']:
                            metrics.append({'locus': locus['locus'], 'repeat': rep+1, 'prediction': kind, **metric})
                statuses.append({'locus': locus['locus'], 'repeat': rep+1, 'state': state, 'n_ALT': len(locus['alts']), 'n_predictors': len(features)})
                if streams and locus['reason'] is None:
                    pred = result['target_compatible'] if result is not None else np.full((len(target_samples),len(locus['alts'])),np.nan)
                    calls = ['./.:' + (','.join(format(v, '.17g') for v in row) if np.isfinite(row).all() else ','.join('.' for _ in row)) for row in pred]
                    streams[rep].write(f"{locus['chrom']}\t{locus['pos']}\t.\t{locus['ref']}\t{','.join(locus['alts'])}\t.\t.\tEND={locus['end']}\tGT:DS" + ''.join('\t'+c for c in calls) + '\n')
    finally:
        for stream in streams:
            stream.close()
    write_tsv(out/'status.tsv', statuses, ['locus','repeat','state','n_ALT','n_predictors'])
    write_tsv(out/'predictors.tsv', chosen_features, ['locus','n_predictors','predictors'])
    write_tsv(out/'accuracy.tsv', metrics, ['locus','repeat','prediction','allele_id','n_pairs','n_truth_available','pearson_r','mse','pearson_status'])
    pooled = defaultdict(list)
    for row in metrics:
        pooled[(row['locus'],row['allele_id'],row['prediction'])].append(row)
    summaries = []
    for (locus,allele,kind), entries in pooled.items():
        values = [r['pearson_r'] for r in entries if r['pearson_r'] is not None]
        summaries.append({'locus':locus,'allele_id':allele,'prediction':kind,
            'repeats_requested':args.repeats,'repeats_with_defined_r':len(values),
            'mean_signed_r':float(np.mean(values)) if values else None,
            'sd_signed_r':float(np.std(values,ddof=1)) if len(values)>1 else None})
    write_tsv(out/'accuracy_summary.tsv',summaries,['locus','allele_id','prediction','repeats_requested','repeats_with_defined_r','mean_signed_r','sd_signed_r'])
    provenance = {'version': '0.1.0', 'options': vars(args), 'sources': {str(p):checksum(p) for p in [args.snps_vcf,args.targets_vcf,args.predictor_pool,args.target_ids,args.splits_json] if p},
                  'reference_samples':len(ref_samples),'target_samples':len(target_samples),'pool_SNPs':len(predictors),
                  'loci':len(chosen_features),'status_counts':dict(Counter(r['state'] for r in statuses)),
                  'accuracy_filter_applied':False,'QC_filters_applied':False,'originals_modified':False,
                  'limitations':['Diploid called GT only. Split records with incompatible ALT counts are not repaired.',
                    'DS constraints establish arithmetic compatibility, not biological correctness.',
                    'Random individual folds may be optimistic for family/population transfer; provide explicit splits for that use.',
                    'Targets with unresolved positional representation remain in input/status ledgers but have no reconstructed output VCF record.']}
    (out/'run.json').write_text(json.dumps(provenance,indent=2))
    print(json.dumps({k:provenance[k] for k in ['reference_samples','target_samples','pool_SNPs','loci','status_counts']},indent=2))
    return provenance


def legacy_main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snps-vcf',required=True,help='Called biallelic SNPs for references and optional targets')
    p.add_argument('--targets-vcf',required=True,help='Called target SNP/SV variants in reference samples')
    p.add_argument('--predictor-pool',help='Exact IDs or CHROM:POS:REF:ALT, one per line; omitted means all supplied SNPs')
    p.add_argument('--target-ids',help='Restrict evaluation to explicitly listed variant IDs')
    p.add_argument('--chromosome',action='append',help='Repeat for multiple target chromosomes')
    p.add_argument('--predictor-mode',choices=['nearest','chromosome','pool'],default='nearest')
    p.add_argument('--neighbors',type=int,default=64,help='Nearest same-chromosome SNPs, within the supplied pool')
    p.add_argument('--outer-folds',type=int,default=5)
    p.add_argument('--inner-folds',type=int,default=3)
    p.add_argument('--repeats',type=int,default=1)
    p.add_argument('--seed',type=int,default=20261121)
    p.add_argument('--alphas',type=float,nargs='+')
    p.add_argument('--splits-json',help='Explicit reference IDs and nested splits for reproducible/group validation')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--check-only',action='store_true')
    p.add_argument('--out',required=True)
    args=p.parse_args(argv)
    if min(args.neighbors,args.repeats,args.threads)<1:
        p.error('Neighbors, repeats and threads must be positive')
    try:
        run(args)
    except (ValueError,KeyError) as exc:
        p.exit(2,str(exc)+'\n')


def main(argv=None):
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["describe"]:
        parser = argparse.ArgumentParser(prog="sv-impute describe", description="Describe the reference VCF and plot its variants.")
        parser.add_argument("--reference", required=True, help="Coordinate-sorted reference VCF")
        parser.add_argument("--chromosome", action="append", help="VCF contig to describe; may be repeated")
        parser.add_argument("--out", required=True, help="New results directory")
        args = parser.parse_args(argv[1:])
        from .benchmark import run as benchmark
        try:
            benchmark(args, describe_only=True)
        except (OSError, ValueError, KeyError) as exc:
            parser.exit(2, str(exc) + "\n")
        return
    if any(arg.split("=")[0] in ("--snps-vcf", "--targets-vcf") for arg in argv):
        return legacy_main(argv)
    parser = argparse.ArgumentParser(description="Predict ALT dosages from local SNPs in a user-defined LD panel.",
                                     epilog="For a reference-only report, use sv-impute describe --help.")
    parser.add_argument("--reference", required=True, help="Mixed, coordinate-sorted VCF with called SNP and SV genotypes")
    panel = parser.add_mutually_exclusive_group(required=True)
    panel.add_argument("--ld-panel-vcf", "--ld", metavar="VCF",
                       help="LD SNP sites, optionally with genotypes for new samples (--ld is an alias)")
    panel.add_argument("--ld-panel-bim", metavar="BIM",
                       help="PLINK variant map defining LD SNPs for validation in the reference samples")
    parser.add_argument("--out", required=True, help="New results directory")
    parser.add_argument("--neighbors", type=int, nargs="+", default=[1, 2, 4, 8, 16, 32, 64, 128],
                        help="Neighborhood sizes to compare (default: 1 2 4 8 16 32 64 128)")
    parser.add_argument("--targets", choices=["sv-indel", "sv", "snp", "all"], default="sv-indel",
                        help="Default: SVs and sequence-resolved indels of any length")
    parser.add_argument("--chromosome", action="append", help="Restrict to this VCF contig; may be repeated")
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--splits-json", help="Use explicit nested splits in reference sample order")
    parser.add_argument("--threads", type=int, default=1, help="Numerical library threads")
    args = parser.parse_args(argv)
    if min(*args.neighbors, args.repeats, args.threads) < 1:
        parser.error("Window sizes, repeats and threads must be positive.")
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[key] = str(args.threads)
    from .benchmark import run as benchmark
    try:
        benchmark(args)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, str(exc) + "\n")
