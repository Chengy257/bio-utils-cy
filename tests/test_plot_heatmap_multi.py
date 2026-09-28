"""Functional tests for bin/plot_heatmap_multi.R (real ComplexHeatmap).

The input contract mirrors what plot_go_heatmap_batch.sh produces:
matrices named <prefix>.<cluster> (the cluster name is filename field 2
after a '.' split) and a one-label-per-line group file aligned to the
matrix columns.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "ComplexHeatmap", "circlize")


def matrix_text(rows, cols, seed=1):
    import random
    rnd = random.Random(seed)
    lines = ["gene\t" + "\t".join(cols)]
    for r in rows:
        lines.append(r + "\t" + "\t".join(
            f"{rnd.uniform(-2, 2):.4f}" for _ in cols))
    return "\n".join(lines) + "\n"


@unittest.skipUnless(RS, "R with getopt+ComplexHeatmap+circlize not available")
class TestPlotHeatmapMulti(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # relative output paths land in the temp dir
        try:
            return subprocess.run([RS, str(BIN / "plot_heatmap_multi.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def _fixture(self, constant_row=False):
        cols = ["S1", "S2", "S3", "S4"]
        rows_a = [f"A{i}" for i in range(1, 9)]
        rows_b = [f"B{i}" for i in range(1, 7)]
        mat_a = self.write("expr.ClusterA.tsv", matrix_text(rows_a, cols, seed=1))
        mat_b_rows = rows_b + (["Bconst"] if constant_row else [])
        mat_b = self.write("expr.ClusterB.tsv",
                           matrix_text(mat_b_rows, cols, seed=2) if not constant_row else
                           matrix_text(rows_b, cols, seed=2) + "Bconst\t5.0\t5.0\t5.0\t5.0\n")
        fl = self.write("files.txt", f"{mat_a}\n{mat_b}\n")
        groups = self.write("groups.txt", "Control\nControl\nTreat\nTreat\n")
        return fl, groups

    def test_rendes_individual_and_combined_pdfs(self):
        fl, groups = self._fixture()
        proc = self._run("-f", fl, "-g", groups, "-o", "hm")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name in ("hm_ClusterA_Single_Heatmap.pdf",
                     "hm_ClusterB_Single_Heatmap.pdf",
                     "hm_Multi_Heatmap.pdf"):
            p = self.tmp / name
            self.assertTrue(p.exists() and p.stat().st_size > 3000, name)
        self.assertIn("Panels plotted: 2", proc.stderr)

    def test_constant_row_dropped_with_note(self):
        # a constant row gives NaN z-scores and is dropped per panel
        fl, groups = self._fixture(constant_row=True)
        proc = self._run("-f", fl, "-g", groups, "-o", "hm")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Retained 6 of 7 rows", proc.stderr)

    def test_group_count_mismatch_fails(self):
        fl, _ = self._fixture()
        groups = self.write("bad_groups.txt", "Control\nTreat\n")
        proc = self._run("-f", fl, "-g", groups, "-o", "hm")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("must match", proc.stderr)

    def test_missing_input_matrix_fails(self):
        _, groups = self._fixture()
        fl = self.write("files.txt", str(self.tmp / "nope.tsv") + "\n")
        proc = self._run("-f", fl, "-g", groups, "-o", "hm")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not found", proc.stderr)

    def test_viridis_color_scheme(self):
        fl, groups = self._fixture()
        proc = self._run("-f", fl, "-g", groups, "-o", "hm", "--color", "viridis")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Color scheme:   viridis", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
