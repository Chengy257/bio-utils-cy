"""Functional tests for bin/plot_heatmap_single.R (real ComplexHeatmap)."""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "ComplexHeatmap", "circlize", "RColorBrewer")


def matrix_text(rows, cols, seed=1, constant_row=None):
    import random
    rnd = random.Random(seed)
    lines = ["gene\t" + "\t".join(cols)]
    for r in rows:
        if constant_row and r == constant_row:
            lines.append(r + "\t" + "\t".join(["5.0"] * len(cols)))
        else:
            lines.append(r + "\t" + "\t".join(
                f"{rnd.uniform(-2, 2):.4f}" for _ in cols))
    return "\n".join(lines) + "\n"


@unittest.skipUnless(RS, "R with getopt+ComplexHeatmap+circlize+RColorBrewer not available")
class TestPlotHeatmapSingle(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # relative outputs land in the temp dir
        try:
            return subprocess.run([RS, str(BIN / "plot_heatmap_single.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def test_renders_pdf_default(self):
        mat = self.write("expr.tsv", matrix_text([f"G{i}" for i in range(1, 11)],
                                                 ["S1", "S2", "S3", "S4"]))
        proc = self._run("-d", mat)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdf = self.tmp / "heatmap.pdf"  # default prefix
        self.assertTrue(pdf.exists() and pdf.stat().st_size > 3000, proc.stderr)
        self.assertIn("Rows (genes):   10 / 10 retained", proc.stderr)

    def test_constant_row_reported_as_rows(self):
        # regression: the notice counted NA cells, not removed rows
        mat = self.write("expr.tsv", matrix_text(
            [f"G{i}" for i in range(1, 10)] + ["CONST"],
            ["S1", "S2", "S3", "S4"], constant_row="CONST"))
        proc = self._run("-d", mat, "-o", "hm")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("caused 1 row(s) to be removed", proc.stderr)
        self.assertIn("9 / 10 retained", proc.stderr)

    def test_crlf_group_labels(self):
        # regression: \r landed in the legend factor levels
        mat = self.write("expr.tsv", matrix_text([f"G{i}" for i in range(1, 11)],
                                                 ["S1", "S2", "S3", "S4"]))
        grp = self.write("groups.txt", "Control\r\nControl\r\nTreat\r\nTreat\r\n")
        proc = self._run("-d", mat, "-g", grp, "-o", "hm")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Groups: Control, Treat", proc.stderr)

    def test_group_count_mismatch_fails(self):
        mat = self.write("expr.tsv", matrix_text([f"G{i}" for i in range(1, 11)],
                                                 ["S1", "S2", "S3", "S4"]))
        grp = self.write("groups.txt", "Control\nTreat\n")
        proc = self._run("-d", mat, "-g", grp, "-o", "hm")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("must match", proc.stderr)

    def test_column_cluster_space_form(self):
        # the "--opt TRUE" space form is the supported one; the "=" form
        # trips a bug inside getopt itself (condition has length > 1)
        mat = self.write("expr.tsv", matrix_text([f"G{i}" for i in range(1, 11)],
                                                 ["S1", "S2", "S3", "S4"]))
        proc = self._run("-d", mat, "-o", "hm", "--column-cluster", "TRUE")
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_show_row_names_false(self):
        mat = self.write("expr.tsv", matrix_text([f"G{i}" for i in range(1, 11)],
                                                 ["S1", "S2", "S3", "S4"]))
        proc = self._run("-d", mat, "-o", "hm", "--show-row-names", "FALSE")
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_missing_data_file(self):
        proc = self._run("-d", self.tmp / "nope.tsv")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("data file not found", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
