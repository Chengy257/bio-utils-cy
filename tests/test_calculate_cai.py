"""Functional tests for bin/calculate_cai.py.

The CAI package is not installed in this environment, so end-to-end
coverage stops at the dependency check; the filtering logic (biopython
only) is exercised directly, including the symmetric-reference-filter
fix.
"""

import importlib
import sys
import unittest

from common import BIN, ScriptTestCase

sys.path.insert(0, str(BIN))
cai_mod = importlib.import_module("calculate_cai")


class CaiFilterTest(ScriptTestCase):

    def test_is_valid_cds(self):
        self.assertTrue(cai_mod.is_valid_cds("ATGGCCTAA"))
        self.assertTrue(cai_mod.is_valid_cds("atggcctaa"))     # case-insensitive
        self.assertFalse(cai_mod.is_valid_cds("ATGGC"))        # not %3
        self.assertFalse(cai_mod.is_valid_cds("ATGNCCAA"))     # N
        self.assertFalse(cai_mod.is_valid_cds(""))             # empty

    def test_reference_filter_is_symmetric_with_targets(self):
        # Regression: the reference filter used to accept lowercase and
        # N-containing sequences that targets would reject.
        ref = (
            ">ok\nATGGCCTAA\n"
            ">lower\natggcctaa\n"     # kept: uppercase-normalized first
            ">withn\nATGNCCAA\n"      # now skipped (kept before)
            ">badlen\nATGGC\n"        # skipped both before and now
        )
        inp = self.write("ref.fa", ref)
        out = self.tmp / "filtered.fa"
        total, kept = cai_mod.filter_reference(str(inp), str(out))
        self.assertEqual((total, kept), (4, 2))
        text = out.read_text()
        self.assertIn(">ok", text)
        self.assertIn(">lower", text)
        self.assertNotIn(">withn", text)

    def test_target_filter(self):
        tgt = ">ok\nATGGCCTAA\n>bad\nATGNCCAA\n"
        inp = self.write("tgt.fa", tgt)
        out = self.tmp / "filtered.fa"
        total, kept = cai_mod.filter_targets(str(inp), str(out))
        self.assertEqual((total, kept), (2, 1))


class CaiE2ETest(ScriptTestCase):

    def test_missing_cai_dependency_clean_exit(self):
        inp = self.write("t.fa", ">a\nATGGCCTAA\n")
        ref = self.write("r.fa", ">r\nATGGCCTAA\n")
        proc = self.run_script("calculate_cai.py", "-i", inp, "-r", ref,
                               "-o", self.tmp / "o.tsv")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("CAI", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_help_works_without_dependency(self):
        proc = self.run_script("calculate_cai.py", "--help")
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_missing_inputs(self):
        proc = self.run_script("calculate_cai.py", "-i", self.tmp / "no.fa",
                               "-r", self.tmp / "no2.fa", "-o", self.tmp / "o.tsv")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
