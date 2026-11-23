import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pysam

from svimpute.cli import main
from svimpute.profile import ReferenceProfile

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class ReferenceProfileTests(unittest.TestCase):
    def test_example_cli_and_lengths(self):
        with TemporaryDirectory() as tmp:
            out = Path(tmp) / "report"
            main(["describe", "--reference", str(EXAMPLES / "allelic_series.vcf"), "--out", str(out)])
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual((summary["samples"], summary["records"], summary["positional_loci"]), (0, 8, 7))
            self.assertEqual((summary["sv_alt_entries"], summary["sv_entries_with_length"]), (9, 8))
            with (out / "alleles_per_locus.tsv").open() as f:
                counts = {(r["scope"], r["allelic_class"]): int(r["loci"]) for r in csv.DictReader(f, delimiter="\t")}
            for scope, expected in [("All loci", [5, 1, 1, 0]), ("Loci with SVs", [4, 1, 1, 0])]:
                self.assertEqual([counts[scope, k] for k in ("Biallelic","Triallelic","4+ alleles","Unresolved")], expected)
            with gzip.open(out / "sv_alleles.tsv.gz", "rt") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            lengths = {r["id"]: (r["length_bp"], r["length_source_or_reason"]) for r in rows}
            self.assertEqual(lengths["inv1"], ("150", "REFERENCE_SPAN"))
            self.assertEqual(lengths["del1"], ("75", "SVLEN"))
            self.assertEqual(lengths["bnd1"], ("", "BREAKEND_OR_TRANSLOCATION"))
            self.assertEqual([r["length_bp"] for r in rows if r["id"] == "ins3"], ["2", "4", "6"])
            for name in ("alleles_per_locus", "sv_types", "sv_lengths"):
                self.assertTrue((out / (name + ".png")).read_bytes().startswith(b"\x89PNG"))

    def test_split_records_placeholders_and_ambiguous_loci(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            vcf = root / "input.vcf"
            header = [s for s in (EXAMPLES / "allelic_series.vcf").read_text().splitlines() if s.startswith("#")]
            rows = [
                ("demo1", 1, "A", "C", "."), ("demo1", 1, "A", "T", "."), ("demo1", 1, "A", "C", "."),
                ("demo1", 2, "A", "<NON_REF>", "."), ("demo1", 3, "A", "G,<NON_REF>", "."),
                ("demo1", 4, "A", "T", "."), ("demo1", 4, "AC", "A", "."),
                ("demo1", 5, "A", "<DEL>", "END=20"), ("demo1", 5, "A", "<DEL>", "END=30"),
                ("demo1", 6, "A", "<DEL>,<DUP>", "SVLEN=-11,22"),
                ("demo2", 1, "A", "G", ".")]
            data = [f"{c}\t{p}\t.\t{ref}\t{alt}\t.\t.\t{info}" for c,p,ref,alt,info in rows]
            vcf.write_text("\n".join(header + data) + "\n")
            profile = ReferenceProfile(root / "profile")
            with pysam.VariantFile(str(vcf)) as src:
                list(profile.scan(src, ["demo1"]))
            self.assertEqual(profile.classes["All loci"], {"Biallelic": 1, "Triallelic": 2, "Unresolved": 3})
            self.assertEqual(profile.lengths["DUP"][22], 1)
            self.assertEqual(profile.lengths["DEL"][11], 1)
            with gzip.open(profile.output / "unresolved_loci.tsv.gz", "rt") as f:
                self.assertEqual({r["reason"] for r in csv.DictReader(f, delimiter="\t")},
                                 {"NO_DISTINCT_EXPLICIT_ALT", "DIFFERENT_REF_AT_ONE_POSITION",
                                  "SYMBOLIC_ALT_WITH_DIFFERENT_SPANS_OR_LENGTHS"})
