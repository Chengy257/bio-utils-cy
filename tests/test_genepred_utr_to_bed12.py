"""Functional tests for bin/genepred_utr_to_bed12.py.

Hand-written genePred rows with coordinates chosen so the expected
BED12 output can be derived by inspection. Key assertions: minus-strand
5'/3' UTR side selection, multi-block UTRs, bin-prefixed 11-column
input, and non-coding row skipping.
"""

import importlib.util
import os
import unittest
from pathlib import Path

from common import BIN, ScriptTestCase

# Standard 10-column genePred rows.
GP_PLUS = (
    "g_plus\tchr1\t+\t100\t900\t200\t700\t3\t100,300,600,\t200,450,900,\n"
)
GP_MINUS = (
    "g_minus\tchr1\t-\t100\t900\t200\t700\t3\t100,300,600,\t200,450,900,\n"
)
GP_SPAN = (
    "g_span\tchr1\t+\t100\t900\t500\t550\t3\t100,400,700,\t300,600,900,\n"
)
GP_NONCODING = (
    "g_nc\tchr1\t+\t100\t900\t100\t100\t2\t100,500,\t200,900,\n"
)
# Bin-prefixed 11-column UCSC knownGene export (bin=585), same gene as GP_PLUS.
GP_BIN = (
    "585\tg_bin\tchr1\t+\t100\t900\t200\t700\t3\t100,300,600,\t200,450,900,\n"
)


def parse_bed(path):
    rows = []
    for line in path.read_text().splitlines():
        rows.append(line.split("\t"))
    return rows


class GenepredUtrToBed12Test(ScriptTestCase):

    def _run(self, gp_text, utr="5UTR", out="utr.bed12", extra_args=()):
        gp = self.write("ann.gp", gp_text)
        out_path = self.tmp / out
        proc = self.run_script(
            "genepred_utr_to_bed12.py", "-i", gp, "-o", out_path, "--utr", utr,
            *extra_args,
        )
        return proc, out_path

    def test_plus_strand_5utr(self):
        proc, out = self._run(GP_PLUS, utr="5UTR")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = parse_bed(out)
        self.assertEqual(len(rows), 1)
        r = rows[0]
        # 5'UTR on plus strand: exon portion left of cdsStart=200 -> (100,200).
        self.assertEqual(r[:3], ["chr1", "100", "200"])
        self.assertEqual(r[3:6], ["g_plus", "0", "+"])
        self.assertEqual(r[9:], ["1", "100,", "0,"])

    def test_plus_strand_3utr(self):
        proc, out = self._run(GP_PLUS, utr="3UTR")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = parse_bed(out)
        self.assertEqual(len(rows), 1)
        # 3'UTR on plus strand: exon portion right of cdsEnd=700 -> (700,900).
        self.assertEqual(rows[0][:3], ["chr1", "700", "900"])
        self.assertEqual(rows[0][9:], ["1", "200,", "0,"])

    def test_minus_strand_5utr_3utr_sides(self):
        proc, out5 = self._run(GP_MINUS, utr="5UTR", out="m5.bed12")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        proc, out3 = self._run(GP_MINUS, utr="3UTR", out="m3.bed12")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # Minus strand: 5'UTR is right of cdsEnd=700, 3'UTR left of cdsStart=200.
        rows5 = parse_bed(out5)
        rows3 = parse_bed(out3)
        self.assertEqual(rows5[0][:3], ["chr1", "700", "900"])
        self.assertEqual(rows5[0][5], "-")
        self.assertEqual(rows3[0][:3], ["chr1", "100", "200"])
        self.assertEqual(rows3[0][5], "-")

    def test_multi_block_utr(self):
        # cds 500-550 splits exon2; 5'UTR spans exon1 + part of exon2,
        # 3'UTR spans the rest of exon2 + exon3.
        proc, out5 = self._run(GP_SPAN, utr="5UTR", out="s5.bed12")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        r = parse_bed(out5)[0]
        self.assertEqual(r[:3], ["chr1", "100", "500"])
        self.assertEqual(r[9:], ["2", "200,100,", "0,300,"])
        proc, out3 = self._run(GP_SPAN, utr="3UTR", out="s3.bed12")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        r = parse_bed(out3)[0]
        self.assertEqual(r[:3], ["chr1", "550", "900"])
        self.assertEqual(r[9:], ["2", "50,200,", "0,150,"])

    def test_bin_prefixed_11col_input(self):
        proc, out = self._run(GP_BIN, utr="5UTR")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = parse_bed(out)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][:3], ["chr1", "100", "200"])
        self.assertEqual(rows[0][3], "g_bin")

    def test_noncoding_row_skipped(self):
        proc, out = self._run(GP_NONCODING, utr="5UTR")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("non-coding", proc.stderr)
        self.assertEqual(parse_bed(out), [])
        # No-UTR output alone (zero rows, no parse errors) is still exit 0.
        proc, out = self._run(GP_NONCODING, utr="5UTR", out="nc3.bed12")
        self.assertEqual(proc.returncode, 0)

    def test_partial_parse_failure_warns_but_exits_0(self):
        bad = "g_short\tchr1\t+\t100\n"
        proc, out = self._run(GP_PLUS + bad, utr="5UTR")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("failed to parse", proc.stderr)
        rows = parse_bed(out)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][3], "g_plus")

    def test_all_rows_unparseable_exits_1(self):
        bad = "g_short\tchr1\t+\t100\n"
        proc, _ = self._run(bad, utr="5UTR")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("all input lines failed to parse", proc.stderr)

    def test_exon_array_mismatch_is_parse_error(self):
        bad = "g_mm\tchr1\t+\t100\t900\t200\t700\t2\t100,300,\t200,\n"
        proc, _ = self._run(bad, utr="5UTR")
        self.assertEqual(proc.returncode, 1)

    def test_both_output_naming(self):
        gp = GP_PLUS + GP_MINUS
        proc, _ = self._run(gp, utr="both", out="utr.bed12")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        f5 = self.tmp / "utr_5UTR.bed12"
        f3 = self.tmp / "utr_3UTR.bed12"
        self.assertTrue(f5.exists())
        self.assertTrue(f3.exists())
        self.assertEqual(len(parse_bed(f5)), 2)
        self.assertEqual(len(parse_bed(f3)), 2)
        # No double suffix like utr_5UTR_5UTR.bed12.
        self.assertFalse((self.tmp / "utr_5UTR_5UTR.bed12").exists())

    def test_both_output_naming_bed_extension(self):
        proc, _ = self._run(GP_PLUS, utr="both", out="utr.bed")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "utr_5UTR.bed").exists())
        self.assertTrue((self.tmp / "utr_3UTR.bed").exists())

    def test_no_suffix_extension_fallback(self):
        proc, _ = self._run(GP_PLUS, utr="both", out="utr")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "utr_5UTR").exists())
        self.assertTrue((self.tmp / "utr_3UTR").exists())

    def test_threads_zero_rejected(self):
        proc, _ = self._run(GP_PLUS, extra_args=("--threads", "0"))
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("--threads must be >= 1", proc.stderr)

    def test_threads_default_reads_buc_threads(self):
        # Regression: --threads was hardcoded to 1 and ignored BUC_THREADS.
        spec = importlib.util.spec_from_file_location(
            "genepred_utr_to_bed12", Path(BIN) / "genepred_utr_to_bed12.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        saved = os.environ.pop("BUC_THREADS", None)
        try:
            os.environ["BUC_THREADS"] = "3"
            parser = mod.build_parser()
            self.assertEqual(
                parser.parse_args(["-i", "a", "-o", "b"]).threads, 3)
            os.environ["BUC_THREADS"] = ""
            parser = mod.build_parser()
            self.assertEqual(
                parser.parse_args(["-i", "a", "-o", "b"]).threads, 1)
        finally:
            if saved is None:
                os.environ.pop("BUC_THREADS", None)
            else:
                os.environ["BUC_THREADS"] = saved


if __name__ == "__main__":
    unittest.main()
