"""Functional tests for bin/mfuzz_soft_cluster.R (real Mfuzz clustering)."""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "Mfuzz", "Biobase")


def matrix_text(rows, cols, seed=1):
    import random
    rnd = random.Random(seed)
    lines = ["gene\t" + "\t".join(cols)]
    for r in rows:
        lines.append(r + "\t" + "\t".join(
            f"{rnd.uniform(0, 10):.4f}" for _ in cols))
    return "\n".join(lines) + "\n"


@unittest.skipUnless(RS, "R with getopt+Mfuzz+Biobase not available")
class TestMfuzzSoftCluster(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # relative outputs land in the temp dir
        try:
            return subprocess.run([RS, str(BIN / "mfuzz_soft_cluster.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def _matrix(self, n_genes=30, cols=("T0", "T6", "T12", "T24"), name="expr.tsv"):
        return self.write(name, matrix_text([f"G{i}" for i in range(1, n_genes + 1)],
                                            cols))

    def test_full_run_outputs(self):
        mat = self._matrix()
        proc = self._run("-i", mat, "-k", "4", "-o", "mf")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdf = self.tmp / "mf_MfuzzPlot.pdf"
        tsv = self.tmp / "mf_Mfuzz_clusterMembership.tsv"
        self.assertTrue(pdf.exists() and pdf.stat().st_size > 2000, proc.stderr)
        self.assertTrue(tsv.exists(), proc.stderr)
        header = tsv.read_text().splitlines()[0].split("\t")
        self.assertEqual(header[:3], ["Gene", "Cluster", "Membership"])
        self.assertIn("Cluster_4", header)
        self.assertIn("Estimated fuzzifier (m):", proc.stderr)
        self.assertIn("Total: 30 genes", proc.stderr)

    def test_no_stray_rplots_pdf(self):
        # regression: filter.std(visu=TRUE) drew before a device existed
        mat = self._matrix()
        proc = self._run("-i", mat, "-k", "3", "-o", "mf")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertFalse((self.tmp / "Rplots.pdf").exists(),
                         "stray Rplots.pdf created")

    def test_reproducible_with_seed(self):
        mat = self._matrix()
        out1 = self._run("-i", mat, "-k", "3", "-o", "r1", "-S", "7")
        out2 = self._run("-i", mat, "-k", "3", "-o", "r2", "-S", "7")
        self.assertEqual(out1.returncode, 0, out1.stderr)
        self.assertEqual(out2.returncode, 0, out2.stderr)
        t1 = (self.tmp / "r1_Mfuzz_clusterMembership.tsv").read_text()
        t2 = (self.tmp / "r2_Mfuzz_clusterMembership.tsv").read_text()
        self.assertEqual(t1, t2)

    def test_duplicate_columns_rejected(self):
        mat = self.write("dup.tsv", matrix_text(
            [f"G{i}" for i in range(1, 11)], ["T0", "T6", "T6", "T24"]))
        proc = self._run("-i", mat, "-k", "3", "-o", "mf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("duplicated sample column names", proc.stderr)

    def test_k_below_two_rejected(self):
        mat = self._matrix()
        proc = self._run("-i", mat, "-k", "1", "-o", "mf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn(">= 2", proc.stderr)

    def test_too_few_genes_after_filtering(self):
        # 3 genes but 4 clusters -> must fail with the filter message
        mat = self._matrix(n_genes=3)
        proc = self._run("-i", mat, "-k", "4", "-o", "mf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("must be >= number of clusters", proc.stderr)

    def test_missing_input(self):
        proc = self._run("-i", self.tmp / "nope.tsv", "-k", "3", "-o", "mf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Input file not found", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
