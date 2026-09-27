"""Functional tests for bin/plot_go_tree_heatmap.R.

The OrgDb needed for real semantic-similarity trees is not installed,
so script-level tests cover the barplot-only paths (which exercise the
description panel and the aplot assembly). Tree-mode assembly was
verified separately with a fake similarity matrix (see phase 7 report).
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "GOSemSim", "ape", "ggtree", "ggplot2", "aplot")

HEADER = "GO_ID\tDescription\tCluster\tNES\tpvalue\n"


@unittest.skipUnless(RS, "R with getopt+GOSemSim+ape+ggtree+ggplot2+aplot not available")
class TestPlotGoTreeHeatmap(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        return subprocess.run([RS, str(BIN / "plot_go_tree_heatmap.R")]
                              + [str(a) for a in args],
                              capture_output=True, text=True, timeout=600, env=env)

    def _tsv(self, rows):
        return self.write("go.tsv", HEADER + rows)

    def test_barplot_only_renders(self):
        # 2 GO terms: no tree, no OrgDb needed; regression for the
        # discrete-y + continuous-ylim crash that made any output
        # impossible, and for the OrgDb hard dependency in this mode.
        tsv = self._tsv("GO:1\tterm one\tC1\t2.5\t0.001\n"
                        "GO:2\tterm two\tC1\t-1.2\t0.040\n")
        out = self.tmp / "out.pdf"
        proc = self._run("-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(out.exists() and out.stat().st_size > 2000, proc.stderr)
        self.assertIn("barplot + description (no tree)", proc.stderr)
        self.assertIn("Plot saved", proc.stderr)

    def test_tree_mode_without_orgdb_degrades(self):
        # >= 3 GO terms wants a tree; without the OrgDb package the
        # script must fall back to barplot-only with a warning, not die.
        tsv = self._tsv("GO:1\tt1\tC1\t2.5\t0.001\n"
                        "GO:2\tt2\tC1\t-1.2\t0.040\n"
                        "GO:3\tt3\tC1\t0.8\t0.010\n")
        out = self.tmp / "out.pdf"
        proc = self._run("-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("falling back to barplot-only", proc.stderr)
        self.assertTrue(out.exists() and out.stat().st_size > 2000, proc.stderr)

    def test_zero_pvalue_clipped_with_warning(self):
        tsv = self._tsv("GO:1\tt1\tC1\t2.5\t0\n"
                        "GO:2\tt2\tC1\t-1.2\t0.040\n")
        out = self.tmp / "out.pdf"
        proc = self._run("-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("clipped", proc.stderr)
        self.assertTrue(out.exists())

    def test_duplicate_go_conflicting_description_warns(self):
        tsv = self._tsv("GO:1\tdesc A\tC1\t2.5\t0.001\n"
                        "GO:1\tdesc B\tC1\t1.5\t0.010\n"
                        "GO:2\tt2\tC1\t-1.2\t0.040\n")
        out = self.tmp / "out.pdf"
        proc = self._run("-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("multiple descriptions", proc.stderr)
        self.assertTrue(out.exists())

    def test_missing_required_columns(self):
        tsv = self.write("bad.tsv", "GO_ID\tNES\nGO:1\t1.0\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "out.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("missing required columns", proc.stderr)

    def test_non_numeric_nes_rejected(self):
        tsv = self._tsv("GO:1\tt1\tC1\thigh\t0.001\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "out.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("non-numeric", proc.stderr)

    def test_missing_input_file(self):
        proc = self._run("-i", self.tmp / "nope.tsv", "-o", self.tmp / "out.pdf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("input file not found", proc.stderr)

    def test_invalid_ontology_rejected(self):
        tsv = self._tsv("GO:1\tt1\tC1\t2.5\t0.001\n")
        proc = self._run("-i", tsv, "-o", self.tmp / "out.pdf", "-O", "XX")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("must be one of BP/MF/CC", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)
        self.assertIn("-log10(p)", proc.stdout)  # documented colour scale


if __name__ == "__main__":
    unittest.main()
