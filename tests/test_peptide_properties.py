"""Functional tests for bin/peptide_properties.py."""

import importlib.util
import os
from pathlib import Path

from common import BIN, ScriptTestCase


class PeptidePropertiesTest(ScriptTestCase):
    FASTA = ">pep1\nPEPTIDEK\n>pep2\nGFLG\n>pep3\nFLPFLK\n"

    def read_table(self, out):
        lines = out.read_text().rstrip("\n").split("\n")
        header = lines[0].split("\t")
        return header, [dict(zip(header, l.split("\t"))) for l in lines[1:]]

    def test_basic_properties(self):
        fa = self.write("p.fa", self.FASTA)
        out = self.tmp / "props.tsv"
        proc = self.run_script("peptide_properties.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        header, rows = self.read_table(out)
        for col in ("id", "length", "molecular_weight", "net_charge_pH7",
                    "isoelectric_point", "average_hydropathy", "AA_A", "AA_K"):
            self.assertIn(col, header)
        self.assertEqual([r["id"] for r in rows], ["pep1", "pep2", "pep3"])
        self.assertEqual(rows[0]["length"], "8")
        # GFLG: ~392 Da (4 residues + water); sanity range
        self.assertGreater(float(rows[1]["molecular_weight"]), 380.0)
        self.assertLess(float(rows[1]["molecular_weight"]), 405.0)

    def test_row_order_deterministic(self):
        # regression: rows were emitted in thread-completion order
        fa = self.write("p.fa", self.FASTA)
        out1 = self.tmp / "run1.tsv"
        out2 = self.tmp / "run2.tsv"
        self.assertEqual(self.run_script("peptide_properties.py", "-i", fa, "-o", out1).returncode, 0)
        self.assertEqual(self.run_script("peptide_properties.py", "-i", fa, "-o", out2).returncode, 0)
        self.assertEqual(out1.read_text(), out2.read_text())

    def test_invalid_sequence_skipped_with_warning(self):
        fa = self.write("p.fa", ">good\nPEPTIDEK\n>bad\nPEPXK\n")
        out = self.tmp / "props.tsv"
        proc = self.run_script("peptide_properties.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("bad", proc.stderr)
        _, rows = self.read_table(out)
        self.assertEqual([r["id"] for r in rows], ["good"])

    def test_all_invalid_exit_1(self):
        fa = self.write("p.fa", ">bad\nXYZ\n")
        out = self.tmp / "props.tsv"
        proc = self.run_script("peptide_properties.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 1)

    def test_visualize_pdf_created(self):
        fa = self.write("p.fa", self.FASTA)
        out = self.tmp / "props.tsv"
        proc = self.run_script("peptide_properties.py", "-i", fa, "-o", out, "--visualize")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdf = self.tmp / "props_properties.pdf"
        self.assertTrue(pdf.exists())
        self.assertGreater(pdf.stat().st_size, 1000)
        self.assertTrue(pdf.read_bytes().startswith(b"%PDF"))

    def test_invalid_threads_clean_error(self):
        fa = self.write("p.fa", self.FASTA)
        out = self.tmp / "props.tsv"
        for bad in ("0", "-2"):
            proc = self.run_script("peptide_properties.py", "-i", fa, "-o", out,
                                   "--threads", bad)
            self.assertEqual(proc.returncode, 2, f"threads={bad}")
            self.assertNotIn("Traceback", proc.stderr, f"threads={bad}")

    def test_missing_input_clean_error(self):
        out = self.tmp / "props.tsv"
        proc = self.run_script("peptide_properties.py", "-i", self.tmp / "no.fa", "-o", out)
        self.assertEqual(proc.returncode, 1)

    def test_threads_default_reads_buc_threads(self):
        # Regression: --threads was hardcoded to 4 and ignored BUC_THREADS.
        spec = importlib.util.spec_from_file_location(
            "peptide_properties", Path(BIN) / "peptide_properties.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        saved = os.environ.pop("BUC_THREADS", None)
        try:
            os.environ["BUC_THREADS"] = "6"
            parser = mod.build_parser()
            self.assertEqual(
                parser.parse_args(["-i", "a", "-o", "b"]).threads, 6)
            os.environ["BUC_THREADS"] = ""
            parser = mod.build_parser()
            self.assertEqual(
                parser.parse_args(["-i", "a", "-o", "b"]).threads, 4)
        finally:
            if saved is None:
                os.environ.pop("BUC_THREADS", None)
            else:
                os.environ["BUC_THREADS"] = saved


if __name__ == "__main__":
    unittest.main()
