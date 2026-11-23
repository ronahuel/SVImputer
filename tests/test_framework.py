import argparse
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pysam
from svimpute.core import _curve, _scale, allele_metrics, project_alt, validate_splits
from svimpute.cli import load_vcf, main, make_splits


class NumericalTests(unittest.TestCase):
    def test_predictions_match_independent_normal_equations(self):
        rng = np.random.default_rng(42)
        x = rng.normal(size=(20, 37)); x[2, 3] = np.nan
        y = rng.uniform(0, 2, size=(20, 2)); fit = np.arange(14)
        alphas = np.array([0.1, 10, 1000.])
        got, _, _, _ = _curve(x, y, fit, x[14:], np.ones(20, dtype=bool), alphas)
        mean = np.nanmean(x[fit], axis=0)
        filled = np.where(np.isfinite(x[fit]), x[fit], mean)
        sd = filled.std(axis=0); sd[sd == 0] = 1
        z = (filled - mean) / sd
        test = (np.where(np.isfinite(x[14:]), x[14:], mean) - mean) / sd
        center = z.mean(axis=0); target_mean = y[fit].mean(axis=0)
        for i, alpha in enumerate(alphas):
            xc = z-center
            beta = np.linalg.solve(xc.T@xc+alpha*np.eye(z.shape[1]), xc.T@(y[fit]-target_mean))
            expected = np.clip((test-center)@beta+target_mean, 0, 2)
            np.testing.assert_allclose(got[i],expected,atol=1e-11,rtol=1e-11)

    def test_scaling_uses_training_only(self):
        _, evaluation, state = _scale(np.array([[0.,np.nan],[2.,np.nan]]),np.array([[100.,1.]]))
        np.testing.assert_array_equal(state['mean'],[1,0])
        self.assertEqual(evaluation[0,0],99)

    def test_multiallelic_projection_and_signed_metrics(self):
        values=project_alt(np.array([[2.,2.,2.],[1.,1.,0.],[np.nan,1.,0.]]))
        np.testing.assert_allclose(values[0],[2/3]*3)
        self.assertTrue((values[:2].sum(axis=1)<=2).all())
        self.assertTrue(np.isnan(values[2]).all())
        metric=allele_metrics(np.array([[0.],[1.],[2.]]),np.array([[2.],[1.],[0.]]),['A'])[0]
        self.assertAlmostEqual(metric['pearson_r'],-1)
        constant=allele_metrics(np.zeros((4,1)),np.zeros((4,1)),['A'])[0]
        self.assertIsNone(constant['pearson_r'])

    def test_nested_splits_reject_leakage(self):
        splits=make_splits(12,3,2,15)
        validate_splits(splits,12,3,2)
        splits[0]['evaluate'].append(splits[0]['fit'][0])
        with self.assertRaises(ValueError):validate_splits(splits,12,3,2)


class VCFTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.samples=[f'animal{i}' for i in range(15)]
        self.refs=list(reversed(self.samples[:12]))
        self.snps=self.root/'snps.vcf';self.targets=self.root/'targets.vcf'
        header='##fileformat=VCFv4.2\n##contig=<ID=chrA,length=1000>\n##contig=<ID=chrB,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Called genotype">\n'
        def gt(n):return ['0/0','0/1','1/1'][n%3]
        with self.snps.open('w')as f:
            f.write(header+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(self.samples)+'\n')
            for j,(chrom,pos)in enumerate([('chrA',10),('chrA',20),('chrB',10)]):
                f.write(f'{chrom}\t{pos}\tp{j}\tA\tC\t.\t.\t.\tGT\t'+'\t'.join(gt((i+j)//(j+1))for i in range(15))+'\n')
        with self.targets.open('w')as f:
            f.write(header+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(self.refs)+'\n')
            calls=['0/0','0/1','1/2','2/2']*3
            f.write('chrA\t50\tsv1\tA\tAT,AG\t.\t.\t.\tGT\t'+'\t'.join(calls)+'\n')
            for name,alt in [('bad1','AC'),('bad2','AA')]:
                f.write(f'chrA\t80\t{name}\tA\t{alt}\t.\t.\t.\tGT\t'+'\t'.join(['1/1']*12)+'\n')
            f.write('chrB\t50\tsvB\tA\tAT\t.\t.\t.\tGT\t'+'\t'.join(gt(i)for i in range(12))+'\n')
        self.pool=self.root/'panel.txt';self.pool.write_text('p0\np2\n')

    def tearDown(self):self.temp.cleanup()

    def arguments(self,dest):
        return ['--snps-vcf',str(self.snps),'--targets-vcf',str(self.targets),
            '--predictor-pool',str(self.pool),'--outer-folds','3','--inner-folds','2',
            '--repeats','2','--out',str(dest)]

    def test_end_to_end_pool_chromosome_and_incompatible_locus(self):
        dest=self.root/'run'
        with contextlib.redirect_stdout(io.StringIO()):main(self.arguments(dest)+['--chromosome','chrA'])
        report=json.loads((dest/'run.json').read_text())
        self.assertEqual(report['pool_SNPs'],2)
        self.assertEqual(report['reference_samples'],12);self.assertEqual(report['target_samples'],3)
        self.assertEqual(report['loci'],2)
        self.assertEqual(report['status_counts']['REFERENCE_EXCEEDS_DIPLOID_COPIES'],2)
        self.assertIn('chrA:10:A:C',(dest/'predictors.tsv').read_text())
        self.assertNotIn('chrB:10:A:C',(dest/'predictors.tsv').read_text())
        for path in dest.glob('imputed_repeat*.vcf'):
            with pysam.VariantFile(str(path))as vcf:
                self.assertEqual(list(vcf.header.samples),self.samples[12:])
                records=list(vcf);self.assertEqual(len(records),1)
                self.assertEqual(records[0].alts,('AT','AG'))
                for sample in vcf.header.samples:
                    dosage=records[0].samples[sample]['DS']
                    self.assertTrue(all(v>=0 for v in dosage))
                    # VCF float parser uses float32; allow its representation error.
                    self.assertLessEqual(sum(dosage),2+2e-7)
        self.assertTrue((dest/'accuracy_summary.tsv').exists())

    def test_unmatched_pool_is_explicit_error(self):
        with self.assertRaisesRegex(ValueError,'absent or ambiguous'):
            load_vcf(self.snps,pool=['does_not_exist'])

    def test_check_only_preserves_source_and_reports_incompatible(self):
        before=self.targets.read_bytes();dest=self.root/'check'
        with contextlib.redirect_stdout(io.StringIO()):main(self.arguments(dest)+['--check-only'])
        self.assertEqual(self.targets.read_bytes(),before)
        self.assertIn('REFERENCE_EXCEEDS_DIPLOID_COPIES',(dest/'status.tsv').read_text())
        self.assertFalse(list(dest.glob('imputed*.vcf')))


if __name__=='__main__':unittest.main()
