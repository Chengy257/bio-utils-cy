"""Functional tests for bin/plot_flower.R.

Rendering geometry (radial petals, upright labels) was verified
visually on 4/5/22-sample renders; these tests cover the CLI contract,
input validation and successful PDF production.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "plotrix")


@unittest.skipUnless(RS, "R with getopt+plotrix not available")
class TestPlotFlower(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        return subprocess.run([RS, str(BIN / "plot_flower.R")]
                              + [str(a) for a in args],
                              capture_output=True, text=True, timeout=600, env=env)

    def _tsv(self, rows, header="sample\tcount\n"):
        return self.write("flower.tsv", header + rows)

    def test_renders_pdf_with_summary(self):
        tsv = self._tsv("S1\t12\nS2\t8\nS3\t25\nS4\t3\n")
        out = self.tmp / "flower.pdf"
        proc = self._run("-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(out.exists() and out.stat().st_size > 1000, proc.stderr)
        self.assertIn("PDF written", proc.stderr)
        self.assertIn("Number of samples: 4", proc.stderr)

    def test_start_angle_zero_renders(self):
        tsv = self._tsv("A\t1\nB\t2\nC\t3\n")
        out = self.tmp / "flower.pdf"
        proc = self._run("-i", tsv, "-o", out, "-s", "0")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(out.exists())

    def test_many_samples_render(self):
        # 22 samples: value labels must stay inside the (short) petals
        tsv = self._tsv("".join(f"S{i}\t{i * 3}\n" for i in range(1, 23)))
        out = self.tmp / "flower22.pdf"
        proc = self._run("-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(out.exists() and out.stat().st_size > 1000)

    def test_negative_count_rejected(self):
        tsv = self._tsv("S1\t12\nS2\t-3\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "o.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Negative counts", proc.stderr)

    def test_non_numeric_count_rejected(self):
        tsv = self._tsv("S1\t12\nS2\tmany\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "o.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("non-numeric", proc.stderr)

    def test_single_sample_rejected(self):
        tsv = self._tsv("S1\t12\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "o.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("At least 2 samples", proc.stderr)

    def test_single_column_rejected(self):
        tsv = self.write("one_col.tsv", "sample\nS1\nS2\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "o.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("at least 2 columns", proc.stderr)

    def test_missing_input(self):
        proc = self._run("-i", self.tmp / "nope.tsv", "-o", self.tmp / "o.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Input file not found", proc.stderr)

    def test_bad_text_cex_rejected(self):
        tsv = self._tsv("S1\t1\nS2\t2\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "o.pdf", "-x", "0")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("positive number", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage", proc.stdout)


if __name__ == "__main__":
    unittest.main()
