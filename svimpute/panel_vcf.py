"""Match LD markers against the mixed reference, then stream target loci."""
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import re
from itertools import groupby

import numpy as np
import pysam

from .cli import locus_records


def kind(record):
    alts = record.alts or ()
    if not alts or any(a in ("*", "<NON_REF>", "<*>") for a in alts):
        return "unresolved"
    if "SVTYPE" in record.info or any(a.startswith(("<", ".")) or a.endswith(".")
                                      or "[" in a or "]" in a for a in alts):
        return "sv"
    if len(record.ref) == 1 and all(len(a) == 1 and a in "ACGT" for a in alts) and record.ref in "ACGT":
        return "snp"
    if any(len(a) != len(record.ref) for a in alts):
        return "indel"
    return "other"


def calls(record, samples):
    alts = record.alts or ()
    counts = np.full((len(samples), len(alts)), np.nan)
    states = set()
    for i, sample in enumerate(samples):
        gt = record.samples[sample].get("GT")
        if not gt or all(a is None for a in gt):
            states.add("MISSING_GT")
        elif len(gt) != 2:
            states.add("NON_DIPLOID_GT_UNSUPPORTED")
        elif None in gt:
            states.add("PARTIAL_GT_UNCOMPUTABLE")
        elif any(a < 0 or a > len(alts) for a in gt):
            states.add("GT_ALLELE_INDEX_OUT_OF_RANGE")
        else:
            for j in range(len(alts)):
                counts[i, j] = gt.count(j + 1)
    return counts, states


@dataclass
class Marker:
    chrom: str
    pos: int
    ref: str
    alt: str
    status: str
    new_calls: object = None
    ref_calls: object = None
    flipped: bool = False
    marker_id: str = "."
    reference_pair: tuple = ()
    panel_chrom: str = ""
    reference_id: str = "."
    position_found: bool = False
    reference_snps: int = 0
    allele_orientation: str = "NOT_APPLICABLE"


def chromosome_resolver(contigs, aliases):
    """Map panel labels to reference labels; exact matches precede numeric aliases."""
    names = set(contigs)
    numeric = defaultdict(list)
    def number(name):
        match = re.fullmatch(r"(?:chr)?([0-9]+)", name, re.IGNORECASE)
        return int(match[1]) if match else None

    for name in names:
        if number(name) is not None:
            numeric[number(name)].append(name)
    overrides = {}
    for spec in aliases or ():
        panel, sep, reference = spec.partition("=")
        if not sep or not panel or not reference or "=" in reference:
            raise ValueError("--chrom-map expects PANEL=REFERENCE.")
        if names and reference not in names:
            raise ValueError(f"Unknown reference contig in --chrom-map: {reference}")
        if panel in overrides and overrides[panel] != reference:
            raise ValueError(f"Conflicting --chrom-map entries for {panel}")
        overrides[panel] = reference

    def resolve(label):
        if label in overrides:
            return overrides[label]
        if label in names:
            return label
        candidates = numeric.get(number(label), ())
        return candidates[0] if len(candidates) == 1 else label
    return resolve


def read_panel_vcf(ld_panel, reference_ids, chromosomes, resolve=str):
    markers = []
    reference_set = set(reference_ids)
    with pysam.VariantFile(str(ld_panel)) as src:
        new_ids = [s for s in src.header.samples if s not in reference_set]
        for r in src:
            chrom = resolve(r.contig)
            if chromosomes and chrom not in chromosomes:
                continue
            m = Marker(chrom, r.pos, r.ref, ",".join(r.alts or ()), "ABSENT_IN_REFERENCE", panel_chrom=r.contig)
            m.marker_id = r.id or "."
            markers.append(m)
            if kind(r) != "snp":
                m.status = "LD_RECORD_REQUIRES_SNP"
                continue
            if len(r.alts) != 1:
                m.new_calls = np.full(len(new_ids), np.nan)
                continue
            values, states = calls(r, new_ids)
            if states & {"NON_DIPLOID_GT_UNSUPPORTED", "GT_ALLELE_INDEX_OUT_OF_RANGE"}:
                raise ValueError(f"Unsupported LD genotype at {r.contig}:{r.pos}")
            m.new_calls = values[:, 0]
    return new_ids, markers


def read_panel_bim(path, chromosomes, resolve=str):
    markers = []
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) not in (5, 6):
                raise ValueError(f"{path}, line {line_number}: expected CHR ID [cM] BP A1 A2.")
            chrom, marker_id = fields[:2]
            panel_chrom, chrom = chrom, resolve(chrom)
            bp, a1, a2 = fields[-3:]
            if chromosomes and chrom not in chromosomes:
                continue
            try:
                pos = int(bp)
            except ValueError:
                raise ValueError(f"{path}, line {line_number}: invalid BP coordinate {bp!r}.") from None
            a1, a2 = a1.upper(), a2.upper()
            status = "ABSENT_IN_REFERENCE"
            if pos < 1:
                status = "BIM_POSITION_UNAVAILABLE"
            markers.append(Marker(chrom, pos, a2, a1, status, new_calls=np.empty(0),
                                  marker_id=marker_id, panel_chrom=panel_chrom))
    return [], markers


def match_panel(reference, ld_panel, chromosomes, panel_format="vcf", profile=None, aliases=None,
                candidate_counts=None):
    with pysam.VariantFile(str(reference)) as src:
        reference_ids = list(src.header.samples)
        resolve = chromosome_resolver(src.header.contigs, aliases)
    if not reference_ids:
        raise ValueError("The mixed reference VCF needs samples with called GT.")
    if panel_format == "bim":
        new_ids, markers = read_panel_bim(ld_panel, chromosomes, resolve)
    elif panel_format == "vcf":
        new_ids, markers = read_panel_vcf(ld_panel, reference_ids, chromosomes, resolve)
    else:
        raise ValueError(f"Unknown LD panel format: {panel_format}")
    by_site = defaultdict(list)
    for m in markers:
        if m.status == "ABSENT_IN_REFERENCE":
            by_site[m.chrom, m.pos].append(m)
    counts_by_chrom = defaultdict(Counter)
    labels_by_pos = defaultdict(set)
    if candidate_counts is not None:
        for (chrom, pos), candidates in by_site.items():
            labels_by_pos[pos].update(m.panel_chrom for m in candidates if m.panel_chrom == chrom)
    last_snp = None
    with pysam.VariantFile(str(reference)) as src:
        stream = profile.scan(src, chromosomes) if profile is not None else src
        for r in stream:
            if chromosomes and r.contig not in chromosomes:
                continue
            category = kind(r)
            counts_by_chrom[r.contig][category] += 1
            if candidate_counts is not None and category == "snp" and last_snp != (r.contig, r.pos):
                for label in labels_by_pos.get(r.pos, ()):
                    if label != r.contig:
                        candidate_counts[label, r.contig] += 1
                last_snp = (r.contig, r.pos)
            candidates = by_site.get((r.contig, r.pos), ())
            for m in candidates:
                m.position_found = True
                if category != "snp":
                    if not m.reference_snps:
                        m.status = "NO_REFERENCE_SNP_AT_POSITION"
                    continue
                m.reference_snps += 1
                if m.reference_snps > 1:
                    m.status = "AMBIGUOUS_REFERENCE_MATCH"
                    continue
                if len(r.alts) != 1:
                    m.status = "MULTIALLELIC_REFERENCE_SNP"
                    continue
                values, states = calls(r, reference_ids)
                if states & {"NON_DIPLOID_GT_UNSUPPORTED", "GT_ALLELE_INDEX_OUT_OF_RANGE"}:
                    m.status = "UNSUPPORTED_REFERENCE_SNP_GT"
                    continue
                m.reference_pair = (r.ref, r.alts[0])
                m.reference_id = r.id or "."
                m.ref_calls = values[:, 0]
                # Alleles only orient new-sample genotypes after the position matches.
                if panel_format == "vcf":
                    if (m.ref, m.alt) == m.reference_pair:
                        m.allele_orientation = "SAME"
                    elif (m.alt, m.ref) == m.reference_pair:
                        m.flipped = True
                        m.allele_orientation = "REVERSED"
                        m.new_calls = 2 - m.new_calls
                    else:
                        m.allele_orientation = "UNRESOLVED_NEW_SAMPLE_DOSAGE"
                        m.new_calls[:] = np.nan
                m.status = "MATCHED"
    usable = defaultdict(list)
    for site, candidates in by_site.items():
        matched = [m for m in candidates if m.status == "MATCHED"]
        if matched:
            first = matched[0]
            conflicts = np.zeros(len(first.new_calls), dtype=bool)
            for m in matched[1:]:
                m.status = "DUPLICATE_LD_POSITION"
                both = np.isfinite(first.new_calls) & np.isfinite(m.new_calls)
                conflicts |= both & (first.new_calls != m.new_calls)
                missing = ~np.isfinite(first.new_calls)
                first.new_calls[missing] = m.new_calls[missing]
            if conflicts.any():
                first.allele_orientation = "DISCORDANT_DUPLICATE_LD_CALLS"
                first.new_calls[conflicts] = np.nan
            usable[site[0]].append(first)
    for group in usable.values():
        group.sort(key=lambda m: m.pos)
    for m in markers:
        if m.status == "ABSENT_IN_REFERENCE" and m.chrom not in counts_by_chrom:
            m.status = "NO_REFERENCE_RECORDS_ON_CONTIG"
    return reference_ids, new_ids, markers, usable, counts_by_chrom


def target_loci(reference, ids, panel_sites, chromosomes, selection):
    accepted = {"sv-indel": {"sv", "indel"}, "sv": {"sv"}, "snp": {"snp"},
                "all": {"sv", "indel", "snp", "other", "unresolved"}}[selection]
    with pysam.VariantFile(str(reference)) as src:
        seen_chrom = set()
        last_chrom, last_pos = None, 0
        for (chrom, pos), group in groupby(src, lambda r: (r.contig, r.pos)):
            if chromosomes and chrom not in chromosomes:
                continue
            if chrom != last_chrom:
                if chrom in seen_chrom:
                    raise ValueError("Reference VCF must be grouped by chromosome.")
                seen_chrom.add(chrom)
                last_chrom, last_pos = chrom, 0
            if pos <= last_pos:
                raise ValueError(f"Reference VCF must be sorted by position on {chrom}.")
            last_pos = pos
            rows = []
            for r in group:
                category = kind(r)
                if category not in accepted:
                    continue
                if category == "snp" and (chrom, pos, r.ref, tuple(r.alts)) in panel_sites:
                    continue
                counts, states = calls(r, ids)
                rows.append(dict(chrom=chrom, pos=pos, ref=r.ref, end=r.stop, id=r.id or ".",
                                 alts=list(r.alts or ()), category=category, counts=counts,
                                 sample_status=states))
            if rows:
                locus = next(locus_records(rows))
                if any(row["category"] == "unresolved" for row in rows):
                    locus["reason"] = "UNRESOLVED_OR_SPANNING_ALLELE"
                yield locus


def neighbors(markers, positions, target_pos, width):
    """Walk out from the insertion point. Ties go to the lower coordinate."""
    right = int(np.searchsorted(positions, target_pos))
    left = right - 1
    if right < len(markers) and positions[right] == target_pos:
        right += 1
    selected = []
    while len(selected) < width and (left >= 0 or right < len(markers)):
        if right == len(markers) or (left >= 0 and target_pos - positions[left] <= positions[right] - target_pos):
            selected.append(markers[left])
            left -= 1
        else:
            selected.append(markers[right])
            right += 1
    return selected
