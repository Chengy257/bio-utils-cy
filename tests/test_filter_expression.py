"""End-to-end tests for bin/filter_expression.py.

Fully offline with hand-computed expectations:
  MADs of the fixture rows: g1=1, g2=0, g3 (all zero -> expression-
  filtered), g4=0, g5=1, g6=NaN.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402

FIXTURE = (
    "gene\tS1\tS2\tS3\tS4\tS5\n"
    "g1\t10\t12\t11\t13\t12\n"    # MAD 1, expressed everywhere
    "g2\t1\t1\t1\t1\t1\n"         # MAD 0, expressed everywhere
    "g3\t0\t0\t0\t0\t0\n"         # expression-filtered
    "g4\t5\t0\t6\t0\t0\n"         # MAD 0, exactly min_samples=2 pass
    "g5\t2\t50\t3\t51\t2\n"       # MAD 1
    "g6\t3\tNA\t4\t5\t6\n"        # NaN -> MAD undefined
)


class TestFilterExpression(ScriptTestCase):

    def _run(self, *args):
        return self.run_script("filter_expression.py", *args)

    def _fixture(self, text=FIXTURE):
        return self.write("expr.tsv", text)

    def _out_genes(self, out):
        return {ln.split("\t")[0] for ln in Path(out).read_text().splitlines()[1:]}

    def test_default_filters_keep_variable_expressed_genes(self):
        inp = self._fixture()
        out = self.tmp / "filtered.tsv"
        proc = self._run("-i", inp, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # MAD >= 75th percentile (1.0): g1 and g5 survive; g6 has NaN MAD
        self.assertEqual(self._out_genes(out), {"g1", "g5"})
        self.assertIn("MAD >= 75th percentile", proc.stderr)
        self.assertIn("1 gene(s) have undefined MAD", proc.stderr)

    def test_percentile_zero_keeps_all_expressed(self):
        inp = self._fixture()
        out = self.tmp / "filtered.tsv"
        proc = self._run("-i", inp, "-o", out, "--mad_percentile", "0")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self._out_genes(out), {"g1", "g2", "g4", "g5"})

    def test_skip_mad_only_expression_filter(self):
        inp = self._fixture()
        out = self.tmp / "filtered.tsv"
        proc = self._run("-i", inp, "-o", out, "--skip-mad")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # g6 passes the expression filter (4 samples >= 1); NaN is a
        # missing measurement, not zero
        self.assertEqual(self._out_genes(out), {"g1", "g2", "g4", "g5", "g6"})

    def test_expression_threshold_and_sample_count(self):
        inp = self._fixture()
        out = self.tmp / "filtered.tsv"
        # min_expr 5 in >= 2 samples: g1 yes, g4 yes (5,6), g5 yes (50,51),
        # g6 yes (5,6 — NaN is a missing measurement, not zero), g2 no
        proc = self._run("-i", inp, "-o", out, "--min_expr", "5",
                         "--min_samples", "2", "--skip-mad")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self._out_genes(out), {"g1", "g4", "g5", "g6"})

    def test_nan_rows_warned_in_expression_filter(self):
        inp = self._fixture()
        out = self.tmp / "filtered.tsv"
        # min_expr 3.5: g1..g5 pass on >= 2 samples except none fail; g6's
        # NaN cannot pass (4 samples >= 3.5, NaN not counted -> still 4 >= 2,
        # so g6 survives the expression filter and is only warned+dropped by MAD)
        proc = self._run("-i", inp, "-o", out, "--min_expr", "3.5",
                         "--mad_percentile", "0")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("contain NaN", proc.stderr)

    def test_help_states_percentile_semantics(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("MAD >= this percentile", proc.stdout)
        self.assertNotIn("Keep top N%", proc.stdout)
        self.assertIn("not scaled by 1.4826", proc.stdout)

    def test_invalid_percentile_rejected(self):
        inp = self._fixture()
        proc = self._run("-i", inp, "-o", self.tmp / "o.tsv", "--mad_percentile", "150")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("--mad_percentile must be within [0, 100]", proc.stderr)

    def test_invalid_min_samples_rejected(self):
        inp = self._fixture()
        proc = self._run("-i", inp, "-o", self.tmp / "o.tsv", "--min_samples", "0")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("--min_samples must be >= 1", proc.stderr)

    def test_no_numeric_columns_fails_clearly(self):
        inp = self.write("expr.tsv", "gene\tsampleA\nfoo\thigh\nbar\tlow\n")
        proc = self._run("-i", inp, "-o", self.tmp / "o.tsv", "--skip-mad")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No numeric columns", proc.stderr)

    def test_all_filtered_out_warns(self):
        inp = self._fixture()
        out = self.tmp / "filtered.tsv"
        proc = self._run("-i", inp, "-o", out, "--min_expr", "100")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("All genes filtered out", proc.stderr)
        self.assertEqual(self._out_genes(out), set())


if __name__ == "__main__":
    unittest.main()
