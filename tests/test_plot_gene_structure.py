"""Functional tests for bin/plot_gene_structure.R.

Runs the real rtracklayer GTF import against small synthetic GTFs.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "ggplot2", "rtracklayer", "dplyr")

GTF_OK = '''chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; exon_id "e1";
chr1\tsrc\tCDS\t120\t180\t.\t+\t.\tgene_id "GENE1"; cds_id "c1";
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "GENE1"; exon_id "e2";
chr1\tsrc\texon\t1000\t1100\t.\t-\t.\tgene_id "GENE2"; exon_id "e3";
chr1\tsrc\texon\t1300\t1400\t.\t-\t.\tgene_id "GENE2"; exon_id "e4";
chr1\tsrc\texon\t1500\t1600\t.\t-\t.\tgene_id "GENE2"; exon_id "e5";
'''

GTF_NO_GENE_ID = '''chr1\tsrc\texon\t100\t200\t.\t+\t.\texon_id "e1";
chr1\tsrc\texon\t300\t400\t.\t+\t.\texon_id "e2";
'''

GTF_MULTICHROM_GENE = '''chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "GENE3"; exon_id "e1";
chr2\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "GENE3"; exon_id "e2";
'''

GTF_DIFF_CHROM_GENES = '''chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; exon_id "e1";
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "GENE1"; exon_id "e2";
chr2\tsrc\texon\t1000\t1100\t.\t-\t.\tgene_id "GENE2"; exon_id "e3";
chr2\tsrc\texon\t1300\t1400\t.\t-\t.\tgene_id "GENE2"; exon_id "e4";
'''

GTF_NO_EXON = '''chr1\tsrc\tCDS\t120\t180\t.\t+\t.\tgene_id "GENE1"; cds_id "c1";
'''


@unittest.skipUnless(RS, "R with getopt+ggplot2+rtracklayer+dplyr not available")
class TestPlotGeneStructure(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # relative output paths land in the temp dir
        try:
            return subprocess.run([RS, str(BIN / "plot_gene_structure.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def test_renders_multigene_pdf(self):
        gtf = self.write("ok.gtf", GTF_OK)
        out = self.tmp / "genes.pdf"
        proc = self._run("-g", gtf, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(out.exists() and out.stat().st_size > 3000, proc.stderr)
        self.assertIn("Plotting 2 gene(s)", proc.stderr)

    def test_missing_gene_id_attribute_fails_cleanly(self):
        # regression: cryptic "arguments imply differing number of rows"
        gtf = self.write("bad.gtf", GTF_NO_GENE_ID)
        proc = self._run("-g", gtf, "-o", self.tmp / "bad.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("no 'gene_id' attribute", proc.stderr)

    def test_gene_filter_with_spaces(self):
        gtf = self.write("ok.gtf", GTF_OK)
        out = self.tmp / "one.pdf"
        proc = self._run("-g", gtf, "-o", out, "-s", "GENE1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Plotting 1 gene(s)", proc.stderr)
        self.assertTrue(out.exists())

    def test_gene_filter_missing_entry_warns(self):
        gtf = self.write("ok.gtf", GTF_OK)
        proc = self._run("-g", gtf, "-o", self.tmp / "w.pdf", "-s", "GENE1,NOTEXIST")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("not found: NOTEXIST", proc.stderr)

    def test_all_filter_genes_missing_fails(self):
        gtf = self.write("ok.gtf", GTF_OK)
        proc = self._run("-g", gtf, "-o", self.tmp / "x.pdf", "-s", "NOTEXIST")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("None of the specified gene IDs", proc.stderr)

    def test_single_gene_spanning_chromosomes_warns(self):
        gtf = self.write("multi.gtf", GTF_MULTICHROM_GENE)
        proc = self._run("-g", gtf, "-o", self.tmp / "m.pdf")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("span multiple chromosomes", proc.stderr)

    def test_genes_on_different_chromosomes_do_not_warn(self):
        # regression: the old check warned whenever two genes sat on
        # different chromosomes, which is harmless
        gtf = self.write("diff.gtf", GTF_DIFF_CHROM_GENES)
        proc = self._run("-g", gtf, "-o", self.tmp / "d.pdf")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertNotIn("span multiple chromosomes", proc.stderr)

    def test_no_exon_features_fails(self):
        gtf = self.write("noexon.gtf", GTF_NO_EXON)
        proc = self._run("-g", gtf, "-o", self.tmp / "n.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No exon features found", proc.stderr)

    def test_missing_gtf_file(self):
        proc = self._run("-g", self.tmp / "nope.gtf", "-o", self.tmp / "o.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("GTF file not found", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
