"""Functional tests for bin/number_duplicates.R."""

import unittest

from common import RSCRIPT, ScriptTestCase


@unittest.skipUnless(RSCRIPT, "no working Rscript found")
class NumberDuplicatesTest(ScriptTestCase):
    def test_duplicates_and_unique_numbered(self):
        inp = self.write("in.tsv", "geneA\ngeneB\ngeneA\ngeneC\ngeneB\n")
        out = self.tmp / "out.tsv"
        proc = self.run_script("number_duplicates.R", "-i", inp, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = out.read_text().rstrip("\n").split("\n")
        self.assertEqual(
            [l.split("\t") for l in lines],
            [["geneA", "geneA_1"], ["geneB", "geneB_1"], ["geneA", "geneA_2"],
             ["geneC", "geneC_1"], ["geneB", "geneB_2"]],
        )

    def test_numeric_values_keep_leading_zeros(self):
        # regression: read.table type conversion turned 001 into 1
        inp = self.write("in.tsv", "001\n002\n001\n")
        out = self.tmp / "out.tsv"
        proc = self.run_script("number_duplicates.R", "-i", inp, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = out.read_text().rstrip("\n").split("\n")
        self.assertEqual(lines[0], "001\t001_1")
        self.assertEqual(lines[2], "001\t001_2")

    def test_hash_and_quotes_verbatim(self):
        # regression: comment.char="#" and quote defaults ate data
        inp = self.write("in.tsv", 'val#ue\n"quoted"\nval#ue\n')
        out = self.tmp / "out.tsv"
        proc = self.run_script("number_duplicates.R", "-i", inp, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = out.read_text().rstrip("\n").split("\n")
        self.assertEqual(lines[0].split("\t")[0], "val#ue")
        self.assertEqual(lines[2], "val#ue\tval#ue_2")

    def test_help_exit_zero(self):
        proc = self.run_script("number_duplicates.R", "-h")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Usage:", proc.stdout)

    def test_no_args_clean_error(self):
        # regression: used to print usage then crash with an R error
        proc = self.run_script("number_duplicates.R")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Error in", proc.stderr)

    def test_missing_input_clean_error(self):
        proc = self.run_script("number_duplicates.R", "-i", self.tmp / "nope.tsv",
                               "-o", self.tmp / "out.tsv")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Error in", proc.stderr)
        self.assertIn("not found", proc.stderr)

    def test_empty_file_clean_error(self):
        inp = self.write("empty.tsv", "")
        out = self.tmp / "out.tsv"
        proc = self.run_script("number_duplicates.R", "-i", inp, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Error in", proc.stderr)

    def test_summary_goes_to_stderr(self):
        inp = self.write("in.tsv", "a\n")
        out = self.tmp / "out.tsv"
        proc = self.run_script("number_duplicates.R", "-i", inp, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Processed 1 rows", proc.stderr)
        self.assertEqual(proc.stdout, "")


if __name__ == "__main__":
    unittest.main()
