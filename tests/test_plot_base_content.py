"""Functional tests for bin/plot_base_content.R.

Hand-computed expectations use the periodic sequence
ACGTACGTAC (x30 = 300 bp): G at positions 3,7 and C at 2,6,10 within
each 10 bp period, i.e. 5 GC (50%) and 4 AT (40%) per period.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "Biostrings", "ggplot2")

FASTA = ">seq one\n" + "ACGTACGTAC" * 30 + "\n"


def read_tsv(path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]
    return rows


@unittest.skipUnless(RS, "R with getopt+Biostrings+ggplot2 not available")
class TestPlotBaseContent(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        # outputs are CWD-relative; run from the temp dir
        old_cwd = os.getcwd()
        os.chdir(self.tmp)
        try:
            return subprocess.run([RS, str(BIN / "plot_base_content.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def _fasta(self, text=FASTA, name="bc.fa"):
        return self.write(name, text)

    def test_even_window_values_and_anchoring(self):
        # w=20 keeps the historical anchoring [i-10, i+9]; the window
        # at position 11 covers 1..20 = exactly one period pair.
        fasta = self._fasta()
        out = self.tmp / "GC_seq_one_window20_statistics_table.tsv"
        proc = self._run("-f", fasta, "-b", "GC", "-w", "20")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(len(rows), 281)          # positions 11..291
        self.assertEqual(rows[0]["index"], "11")
        self.assertEqual(rows[0]["Content"], "50")  # 10 GC / 20 bp

    def test_odd_window_spans_full_w(self):
        # regression: odd windows analysed w-1 bp. w=15 at position 8
        # now covers 1..15 = 7 GC / 15 = 46.67%.
        fasta = self._fasta()
        out = self.tmp / "GC_seq_one_window15_statistics_table.tsv"
        proc = self._run("-f", fasta, "-b", "GC", "-w", "15")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(len(rows), 286)            # positions 8..293
        self.assertEqual(rows[0]["index"], "8")
        self.assertEqual(rows[0]["Content"], "46.67")

    def test_at_content(self):
        fasta = self._fasta()
        out = self.tmp / "AT_seq_one_window10_statistics_table.tsv"
        proc = self._run("-f", fasta, "-b", "AT", "-w", "10")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(rows[0]["Content"], "50")  # 5 AT / 10 bp

    def test_cutoff_accepted_and_plot_written(self):
        # regression: "-c 50" used to die in getopt (argflag 0)
        fasta = self._fasta()
        proc = self._run("-f", fasta, "-b", "GC", "-w", "20", "-c", "50")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdf = self.tmp / "GC_seq_one_window20_plot.pdf"
        self.assertTrue(pdf.exists() and pdf.stat().st_size > 1000, proc.stderr)

    def test_output_prefix_multi_sequence_distinct(self):
        # regression: with -o and several sequences every sequence
        # wrote to the same files (only the last survived)
        fasta = self.write("two.fa",
                           ">seqA\n" + "ACGT" * 20 + "\n>seqB\n" + "GGCC" * 20 + "\n")
        proc = self._run("-f", fasta, "-b", "GC", "-w", "16", "-o", "myout")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for seq in ("seqA", "seqB"):
            tsv = self.tmp / f"myout_{seq}_statistics_table.tsv"
            pdf = self.tmp / f"myout_{seq}_plot.pdf"
            self.assertTrue(tsv.exists(), proc.stderr)
            self.assertTrue(pdf.exists(), proc.stderr)
        self.assertEqual(read_tsv(self.tmp / "myout_seqB_statistics_table.tsv")[0]["Content"],
                         "100")  # GGCC... is all GC

    def test_output_prefix_single_sequence_exact(self):
        fasta = self._fasta()
        proc = self._run("-f", fasta, "-b", "GC", "-w", "20", "-o", "solo")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "solo_statistics_table.tsv").exists())
        self.assertTrue((self.tmp / "solo_plot.pdf").exists())

    def test_statistics_table_is_tsv(self):
        # the table was tab-delimited but carried an .xls extension
        fasta = self._fasta()
        proc = self._run("-f", fasta, "-b", "GC", "-w", "20")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertFalse((self.tmp / "GC_seq_one_window20_statistics_table.xls").exists())
        self.assertTrue((self.tmp / "GC_seq_one_window20_statistics_table.tsv").exists())

    def test_window_larger_than_sequence_skipped(self):
        fasta = self.write("short.fa", ">tiny\nACGTACGTAA\n")
        proc = self._run("-f", fasta, "-b", "GC", "-w", "50")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("exceeds sequence length", proc.stderr)

    def test_invalid_base_rejected(self):
        fasta = self._fasta()
        proc = self._run("-f", fasta, "-b", "X", "-w", "20")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Invalid base type", proc.stderr)

    def test_bare_cutoff_rejected(self):
        # regression: a bare "--cutoff" was silently read as TRUE and
        # drew the reference line at 1; getopt now demands a value
        fasta = self._fasta()
        proc = self._run("-f", fasta, "-b", "GC", "-w", "20", "--cutoff")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertNotIn("[INFO] Done", proc.stderr)

    def test_bad_plot_dimensions_rejected(self):
        fasta = self._fasta()
        for flag in ("--plot-width", "--plot-height"):
            proc = self._run("-f", fasta, "-b", "GC", "-w", "20", flag, "0")
            self.assertEqual(proc.returncode, 1, flag)
            self.assertIn("positive number", proc.stderr)

    def test_missing_required_args(self):
        proc = self._run("-b", "GC", "-w", "20")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Missing required arguments", proc.stderr)

    def test_missing_fasta_file(self):
        proc = self._run("-f", self.tmp / "nope.fa", "-b", "GC", "-w", "20")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not found", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage", proc.stdout)


if __name__ == "__main__":
    unittest.main()
