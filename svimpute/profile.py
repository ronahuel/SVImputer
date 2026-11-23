"""Describe the alleles represented in a coordinate-sorted reference VCF."""
from collections import Counter, defaultdict
from contextlib import ExitStack
import gzip
from itertools import groupby
import json
from pathlib import Path

from .report import writer

PLACEHOLDERS = {"*", "<NON_REF>", "<*>", "."}
BINS = ("Biallelic", "Triallelic", "4+ alleles", "Unresolved")
BREAKENDS = {"BND", "TRA", "CTX"}


def info_value(record, name, index, shared=False):
    if name not in record.info:
        return None
    value = record.info.get(name)
    if isinstance(value, tuple):
        if len(value) == len(record.alts or ()):
            return value[index]
        return value[0] if shared and len(value) == 1 else None
    if shared or len(record.alts or ()) == 1:
        return value
    return None


def sv_type(record, alt, index):
    if "[" in alt or "]" in alt or alt.startswith(".") or alt.endswith("."):
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].split(":")[0].upper()
    annotation = info_value(record, "SVTYPE", index, shared=True)
    return str(annotation).upper() if annotation and annotation != "." else None


def sv_length(record, alt, index, svtype):
    if svtype in BREAKENDS:
        return None, "BREAKEND_OR_TRANSLOCATION"
    size = info_value(record, "SVLEN", index)
    if size is not None:
        try:
            length = abs(int(size))
        except (TypeError, ValueError):
            return None, "INVALID_SVLEN"
        if length:
            return length, "SVLEN"
    if not alt.startswith("<") and svtype in {"INS", "DEL"}:
        difference = abs(len(alt) - len(record.ref))
        if difference:
            return difference, "SEQUENCE_LENGTH_DIFFERENCE"
    # In older VCFs an anchored symbolic record can give an affected interval.
    # A shared END cannot resolve individual lengths at a multiallelic site.
    if (alt.startswith("<") and len(record.alts) == 1 and svtype in {"DEL", "DUP", "INV", "CNV"}
            and record.stop > record.pos and record.stop != record.start + len(record.ref)):
        return record.stop - record.pos, "REFERENCE_SPAN"
    return None, "LENGTH_UNAVAILABLE"


class ReferenceProfile:
    def __init__(self, output):
        self.output = Path(output)
        self.chromosomes = defaultdict(Counter)
        self.classes = defaultdict(Counter)
        self.types = Counter()
        self.lengths = defaultdict(Counter)

    def scan(self, vcf, chromosomes=None):
        """Yield the same in-scope records to the panel matcher during its first pass."""
        from .panel_vcf import kind

        self.output.mkdir(parents=True, exist_ok=True)
        self.samples = len(vcf.header.samples)
        self.version = vcf.header.version
        self.requested = list(chromosomes or [])
        scoped = (r for r in vcf if not chromosomes or r.contig in chromosomes)
        seen_contigs = set()
        last_chrom, last_pos = None, 0
        with ExitStack() as stack:
            def gz_table(name, fields):
                handle = stack.enter_context(gzip.open(self.output / name, "wt", newline=""))
                return writer(handle, fields)
            allele_table = gz_table("sv_alleles.tsv.gz",
                ["chrom", "pos", "id", "alt_index", "sv_type", "length_bp", "length_source_or_reason"])
            unresolved = gz_table("unresolved_loci.tsv.gz", ["chrom", "pos", "records", "reason"])
            for (chrom, pos), records in groupby(scoped, lambda r: (r.contig, r.pos)):
                if last_chrom != chrom:
                    if chrom in seen_contigs:
                        raise ValueError(f"Reference VCF is not grouped by chromosome at {chrom}:{pos}")
                    seen_contigs.add(chrom)
                    last_chrom, last_pos = chrom, 0
                if pos <= last_pos:
                    raise ValueError(f"Reference positions are out of order at {chrom}:{pos}")
                last_pos = pos
                refs, alts, symbols = set(), set(), {}
                reason, has_sv, rows = "", False, 0
                counts = self.chromosomes[chrom]
                for record in records:
                    rows += 1
                    counts[kind(record)] += 1
                    refs.add(record.ref)
                    for i, alt in enumerate(record.alts or ()):
                        if alt in PLACEHOLDERS:
                            continue
                        alts.add(alt)
                        svtype = sv_type(record, alt, i)
                        if not svtype:
                            continue
                        has_sv = True
                        self.types[svtype] += 1
                        length, source = sv_length(record, alt, i, svtype)
                        if length is not None:
                            self.lengths[svtype][length] += 1
                        allele_table.writerow([chrom, pos, record.id or ".", i+1, svtype, length, source])
                        if alt.startswith("<"):
                            identity = (record.stop, length)
                            if alt in symbols and symbols[alt] != identity:
                                reason = "SYMBOLIC_ALT_WITH_DIFFERENT_SPANS_OR_LENGTHS"
                            symbols[alt] = identity
                    yield record
                alts.difference_update(refs)
                if len(refs) != 1:
                    reason = "DIFFERENT_REF_AT_ONE_POSITION"
                elif not alts:
                    reason = "NO_DISTINCT_EXPLICIT_ALT"
                if reason:
                    group = "Unresolved"
                    unresolved.writerow([chrom, pos, rows, reason])
                else:
                    group = BINS[min(len(alts), 3)-1]
                counts["loci"] += 1
                self.classes["All loci"][group] += 1
                if has_sv:
                    self.classes["Loci with SVs"][group] += 1

    def save(self):
        from .profile_plots import draw_figures
        total = sum(self.chromosomes.values(), Counter())
        categories = ("snp", "sv", "indel", "other", "unresolved")
        by_type = []
        for name, count in sorted(self.types.items()):
            known = sum(self.lengths[name].values())
            by_type.append((name, count, known, count-known))
        with (self.output / "sv_types.tsv").open("w", newline="") as f:
            writer(f, ["type", "alt_entries", "known_length", "unknown_length"]).writerows(by_type)
        with (self.output / "alleles_per_locus.tsv").open("w", newline="") as f:
            table = writer(f, ["scope", "allelic_class", "loci"])
            for scope in ("All loci", "Loci with SVs"):
                table.writerows((scope, label, self.classes[scope][label]) for label in BINS)
        with (self.output / "chromosomes.tsv").open("w", newline="") as f:
            table = writer(f, ["chrom", "loci", *categories])
            for chrom, counts in self.chromosomes.items():
                table.writerow([chrom, counts["loci"], *(counts[c] for c in categories)])
        records = {name: total[name] for name in categories}
        stats = dict(vcf_version=self.version, samples=self.samples, selected_chromosomes=self.requested,
                     chromosomes_with_records=len(self.chromosomes), records=sum(records.values()),
                     records_by_category=records, positional_loci=total["loci"],
                     loci_with_svs=sum(self.classes["Loci with SVs"].values()),
                     sv_alt_entries=sum(self.types.values()), sv_entries_with_length=sum(r[2] for r in by_type),
                     sv_entries_without_length=sum(r[3] for r in by_type))
        (self.output / "summary.json").write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
        draw_figures(self.output, self.classes, self.types, self.lengths)
        print(f"Reference VCF: {self.samples} samples, {stats['records']} records, "
              f"{stats['positional_loci']} loci, {stats['loci_with_svs']} loci containing SVs.")
        return stats
