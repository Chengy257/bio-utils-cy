"""Functional tests for bin/bed12_effective_length.py.

Pure-Python interval merging replaces the former pybedtools dependency,
so these tests run without external tools.
"""

import unittest

from common import ScriptTestCase


class Bed12EffectiveLengthTest(ScriptTestCase):

    def _run(self, bed_text, extra_args=()):
        bed = self.write("in.bed", bed_text)
        proc = self.run_script("bed12_effective_length.py", "-i", bed, *extra_args)
        return proc

    def test_overlapping_merge(self):
        proc = self._run("chr1\t100\t200\t.\t0\t+\nchr1\t150\t300\t.\t0\t+\n")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Effective length: 200 bp", proc.stdout)

    def test_bookended_and_disjoint(self):
        # Bookended [100,200)+[200,300) merge to 200; disjoint [500,600) adds 100.
        proc = self._run(
            "chr1\t100\t200\ta\t0\t+\nchr1\t200\t300\tb\t0\t+\nchr1\t500\t600\tc\t0\t+\n"
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Effective length: 300 bp", proc.stdout)

    def test_chromosomes_merged_independently(self):
        proc = self._run("chr1\t100\t200\ta\nchr2\t150\t300\tb\n")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Effective length: 250 bp", proc.stdout)

    def test_unsorted_input(self):
        # Input not sorted: merge must still be complete (former bedtools
        # merge required pre-sorted input).
        proc = self._run(
            "chr1\t500\t600\tc\nchr1\t100\t250\ta\nchr1\t150\t300\tb\n"
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Effective length: 300 bp", proc.stdout)  # 100-300 + 500-600 = 200+100

    def test_bed12_full_output_columns(self):
        bed = (
            "chr1\t100\t200\tname1\t0\t+\t100\t200\t0\t2\t50,50,\t0,50,\n"
            "chr1\t150\t300\tname2\t0\t+\t150\t300\t0\t1\t150,\t0,\n"
        )
        out = self.tmp / "merged.bed"
        proc = self._run(bed, ["-o", out])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Effective length: 200 bp", proc.stdout)
        rows = out.read_text().splitlines()
        self.assertEqual(len(rows), 1)
        # Remaining columns come from the first interval of the group.
        self.assertEqual(
            rows[0].split("\t"),
            ["chr1", "100", "300", "name1", "0", "+", "100", "200", "0", "2", "50,50,", "0,50,"],
        )

    def test_malformed_lines_skipped(self):
        bed = (
            "chr1\t100\t200\tok\n"
            "chr1\tXX\t300\tbad\n"
            "chr1\t400\n"          # too few columns
            "# comment\n"
            "track name=x\n"
            "chr1\t300\t400\talso_ok\n"
        )
        proc = self._run(bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("malformed", proc.stderr)
        # The bad line is skipped, so only [100,200) and [300,400) remain.
        self.assertIn("Effective length: 200 bp", proc.stdout)

    def test_empty_and_invalid_input_exit_1(self):
        proc = self._run("# only a comment\n")
        self.assertEqual(proc.returncode, 1)

        proc = self._run("chr1\tXX\t300\tbad\n")
        self.assertEqual(proc.returncode, 1)

    def test_missing_input_exit_1(self):
        proc = self.run_script("bed12_effective_length.py", "-i", self.tmp / "nope.bed")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
