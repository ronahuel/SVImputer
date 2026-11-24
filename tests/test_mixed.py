import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from svimpute.cli import main
from svimpute.report import Accuracy

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class MixedInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def run_cli(self, name, ld=None, extra=(), reference=None, panel_option="--ld"):
        dest = self.root / name
        main(["--reference", str(reference or EXAMPLES / "reference.vcf"),
              panel_option, str(ld or EXAMPLES / "ld.vcf"), "--out", str(dest),
              "--outer-folds", "3", "--inner-folds", "2", "--neighbors", "1", "2", "4", *extra])
        return dest

    def test_mixed_reference_chromosomes_predictions_and_peak(self):
        dest = self.run_cli("run")
        counts = {r["chrom"]: r for r in records(dest / "chromosomes.tsv")}
        self.assertEqual([(counts[c]["reference_snps"], counts[c]["ld_matched"]) for c in "ABC"],
                         [("5", "4"), ("1", "1"), ("0", "0")])
        rows = records(dest / "accuracy.tsv.gz")
        self.assertEqual({r["reason"] for r in rows if r["chrom"] == "C"}, {"NO_LOCAL_LD_SNPS"})
        self.assertTrue(all(r["used"] == "1" for r in rows if r["chrom"] == "B"))
        windows = [r for r in records(dest / "window_accuracy.tsv") if r["scope"] == "genome"]
        peak = max(float(r["mean_r_paired"]) for r in windows)
        expected = [int(r["k"]) for r in windows if abs(float(r["mean_r_paired"]) - peak) < 1e-12]
        self.assertEqual(json.loads((dest / "run.json").read_text())["best_caps"], expected)
        from matplotlib.image import imread
        self.assertGreater(imread(dest / "accuracy_by_neighbors.png").shape[0], 100)
        new = records(dest / "dosage_ld.tsv.gz")
        self.assertTrue(new and all(0 <= float(r["new01"]) <= 2 for r in new))
        self.assertEqual(next(r for r in records(dest / "ld_matches.tsv") if r["pos"] == "20")["reversed_alleles"], "True")

    def test_sites_only_gzip_computes_the_same_heldout_dosages(self):
        inp = (EXAMPLES / "ld.vcf").read_text().splitlines()
        sites = self.root / "sites.vcf.gz"
        import pysam
        with pysam.BGZFile(str(sites), "w") as f:
            f.write(("\n".join(line if line.startswith("##") else "\t".join(line.split("\t")[:8])
                               for line in inp) + "\n").encode())
        base = self.run_cli("genotyped", extra=["--chromosome", "A"])
        other = self.run_cli("positions", sites, ["--chromosome", "A"])
        self.assertEqual(records(base / "dosage_cv.tsv.gz"), records(other / "dosage_cv.tsv.gz"))
        self.assertFalse((other / "dosage_ld.tsv.gz").exists())

    def test_ambiguous_and_mismatched_ld_and_untyped_target(self):
        ld = self.root / "ld.vcf"
        lines = (EXAMPLES / "ld.vcf").read_text().splitlines()
        one = next(line for line in lines if line.startswith("A\t10\t"))
        lines += [one, one.replace("A\t10\t", "A\t40\t").replace("\tA\tC\t", "\tA\tT\t"),
                  one.replace("A\t10\t", "A\t90\t")]
        ld.write_text("\n".join(lines) + "\n")
        dest = self.run_cli("other", ld, ["--chromosome", "A", "--targets", "snp"])
        statuses = [r["status"] for r in records(dest / "ld_matches.tsv")]
        self.assertEqual(statuses.count("DUPLICATE_LD_POSITION"), 2)
        self.assertIn("UNRESOLVED_NEW_SAMPLE_DOSAGE",
                      {r["allele_orientation"] for r in records(dest / "ld_matches.tsv")})
        self.assertIn("ABSENT_IN_REFERENCE", statuses)
        for feature in records(dest / "neighbors.tsv.gz"):
            self.assertNotIn(feature["pos"], feature["ld_positions_nearest_first"].split(","))
        self.assertIn("SNP_A", {r["id"] for r in records(dest / "accuracy.tsv.gz")})

    def test_unknown_chromosome_never_gets_a_best_window(self):
        dest = self.run_cli("absent", extra=["--chromosome", "unknown"])
        self.assertEqual(json.loads((dest / "run.json").read_text())["best_caps"], [])
        self.assertEqual(records(dest / "chromosomes.tsv")[0]["ld_matched"], "0")

    def test_multiallelic_and_missing_predictors(self):
        src = (EXAMPLES / "reference.vcf").read_text().splitlines()
        target = next(i for i, v in enumerate(src) if "\tINS_A\t" in v)
        fields = src[target].split("\t")
        fields[4] = "AT,AG"
        fields[9:] = ["0/0", "1/2", "2/2"] * 12
        src[target] = "\t".join(fields)
        reference = self.root / "multi.vcf"
        reference.write_text("\n".join(src) + "\n")
        ld = self.root / "masked.vcf"
        ld.write_text("\n".join(v if v.startswith("#") else v.rsplit("\t", 1)[0] + "\t./."
                                for v in (EXAMPLES / "ld.vcf").read_text().splitlines()) + "\n")
        dest = self.run_cli("multi", ld, ["--chromosome", "A"], reference)
        rows = records(dest / "dosage_ld.tsv.gz")
        self.assertEqual({r["alt"] for r in rows}, {"AT", "AG"})
        for k in {"1", "2", "4"}:
            self.assertLessEqual(sum(float(r["new01"]) for r in rows if r["k"] == k), 2.0000001)
        self.assertTrue(all(r["new02"] == "." for r in rows))

    def test_same_alleles_compared_across_windows(self):
        acc = Accuracy()
        acc.add_locus("A", {1: {(1, 0): -.5, (1, 1): .9}, 4: {(1, 0): -.2}}, {1: 1, 4: 3})
        rows = [r for r in acc.rows() if r["scope"] == "genome"]
        np.testing.assert_allclose([r["mean_r_available"] for r in rows], [.2, -.2])
        np.testing.assert_allclose([r["mean_r_paired"] for r in rows], [-.5, -.2])
        self.assertEqual([r["k"] for r in rows if r["best"]], [4])

    def test_bim_and_vcf_have_identical_reference_predictions(self):
        vcf = self.run_cli("vcf", panel_option="--ld-panel-vcf")
        bim = self.run_cli("bim", EXAMPLES / "ld.bim", panel_option="--ld-panel-bim")
        for filename in ("chromosomes.tsv", "neighbors.tsv.gz", "dosage_cv.tsv.gz",
                         "accuracy.tsv.gz", "window_accuracy.tsv"):
            self.assertEqual(records(vcf / filename), records(bim / filename), filename)
        self.assertEqual(json.loads((bim / "run.json").read_text())["ld_panel_format"], "bim")
        self.assertEqual({row["status"] for row in records(bim / "ld_matches.tsv")}, {"MATCHED"})
        self.assertFalse((bim / "dosage_ld.tsv.gz").exists())
        self.assertTrue((vcf / "dosage_ld.tsv.gz").exists())

    def test_bim_missing_alleles_and_coordinate_ambiguity(self):
        from svimpute.panel_vcf import match_panel
        lines = (EXAMPLES / "reference.vcf").read_text().splitlines()
        i = next(i for i, line in enumerate(lines) if line.startswith("A\t10\t"))
        lines.insert(i+1, lines[i].replace("\tA\tC\t", "\tA\tG\t"))
        reference = self.root / "mixed.vcf"
        reference.write_text("\n".join(lines) + "\n")
        bim = self.root / "map.bim.gz"
        with gzip.open(bim, "wt") as stream:
            stream.write("A ambiguous 0 10 0 .\nA known 20 0 c\nA absent 0 999 A T\n"
                         "A wrong 0 30 A T\nA unknown_bp 0 0 A C\nA repeated_bp 0 40 A C\n"
                         "A repeated_bp2 0 40 C A\nB unambiguous 0 50 0 .\n")
        _, new, markers, usable, _ = match_panel(reference, bim, None, "bim")
        self.assertEqual(new, [])
        self.assertEqual([m.status for m in markers],
                         ["AMBIGUOUS_REFERENCE_MATCH", "MATCHED", "ABSENT_IN_REFERENCE",
                          "MATCHED", "BIM_POSITION_UNAVAILABLE", "MATCHED",
                          "DUPLICATE_LD_POSITION", "MATCHED"])
        self.assertEqual(usable["A"][0].reference_pair, ("A", "C"))

    def test_position_only_pool_ignores_bim_ids_and_both_alleles(self):
        inp = self.root / "panel.bim"
        inp.write_text("7 unused 0 10 T G\n7 other 0 20 0 0\n7 absent 0 999 T A\n")
        out = self.run_cli("audit", inp, ["--chrom-map", "7=A", "--chromosome", "A", "--overlap-only"],
                           panel_option="--ld-panel-bim")
        result = json.loads((out/"run.json").read_text())
        self.assertEqual((result["shared_snp_positions"], result["matched_ld_snps"]), (2, 2))
        self.assertFalse((out / "accuracy.tsv.gz").exists())
        rows = records(out / "ld_matches.tsv")
        self.assertEqual([(r["panel_chrom"], r["chrom"], r["status"]) for r in rows],
                         [("7", "A", "MATCHED"), ("7", "A", "MATCHED"), ("7", "A", "ABSENT_IN_REFERENCE")])
        pending = self.run_cli("unmapped", inp, ["--overlap-only"], panel_option="--ld-panel-bim")
        self.assertEqual(records(pending/"chromosome_candidates.tsv")[0],
                         {"panel_chrom":"7", "candidate_reference_chrom":"A", "coincident_snp_positions":"2"})
        self.assertEqual(json.loads((pending/"run.json").read_text())["matched_ld_snps"], 0)

    def test_automatic_numeric_alias_requires_a_unique_reference_name(self):
        from svimpute.panel_vcf import chromosome_resolver
        self.assertEqual(chromosome_resolver(["chr7"], None)("07"), "chr7")
        self.assertEqual(chromosome_resolver(["chr7", "7"], None)("07"), "07")
        self.assertEqual(chromosome_resolver(["chr7", "7"], None)("7"), "7")

    def test_bim_and_vcf_are_mutually_exclusive(self):
        from contextlib import redirect_stderr
        from io import StringIO
        for extra in ([], ["--ld-panel-vcf", "a.vcf", "--ld-panel-bim", "b.bim"],
                      ["--ld", "a.vcf", "--ld-panel-bim", "b.bim"]):
            with self.subTest(extra=extra), redirect_stderr(StringIO()), self.assertRaises(SystemExit) as error:
                main(["--reference", "reference.vcf", "--out", str(self.root / "no_output"), *extra])
            self.assertEqual(error.exception.code, 2)
        self.assertFalse((self.root / "no_output").exists())
