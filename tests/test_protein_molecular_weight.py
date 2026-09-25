"""Functional tests for bin/protein_molecular_weight.py."""

import csv

from common import ScriptTestCase


class ProteinMolecularWeightTest(ScriptTestCase):
    def read_rows(self, out):
        return {r["ID"]: r for r in csv.DictReader(out.open(), delimiter="\t")}

    def test_known_monoisotopic_weight(self):
        # GFLG: 57.02146+147.06841+113.08406+57.02146 + 18.01056 = 392.20595
        fa = self.write("p.fa", ">pep\nGFLG\n")
        out = self.tmp / "w.tsv"
        proc = self.run_script("protein_molecular_weight.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = self.read_rows(out)
        self.assertEqual(rows["pep"]["Weight_Da"], "392.21")
        self.assertEqual(rows["pep"]["Weight_kDa"], "0.39")
        self.assertEqual(rows["pep"]["Length"], "4")

    def test_lowercase_normalized(self):
        fa = self.write("p.fa", ">pep\ngflg\n")
        out = self.tmp / "w.tsv"
        proc = self.run_script("protein_molecular_weight.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)["pep"]["Weight_Da"], "392.21")

    def test_all_invalid_sequence_reports_na(self):
        # regression: an all-invalid sequence used to report ~18.01 Da
        fa = self.write("p.fa", ">bad\nXXX\n")
        out = self.tmp / "w.tsv"
        proc = self.run_script("protein_molecular_weight.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("no standard residues", proc.stderr)
        rows = self.read_rows(out)
        self.assertEqual(rows["bad"]["Weight_Da"], "NA")
        self.assertEqual(rows["bad"]["Length"], "0")

    def test_partial_invalid_warns(self):
        fa = self.write("p.fa", ">mix\nAXA\n")
        out = self.tmp / "w.tsv"
        proc = self.run_script("protein_molecular_weight.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("1 non-standard residue", proc.stderr)
        rows = self.read_rows(out)
        # 71.03711*2 + 18.01056 = 160.08478
        self.assertEqual(rows["mix"]["Weight_Da"], "160.08")
        self.assertEqual(rows["mix"]["Length"], "2")

    def test_sequence_column_verbatim(self):
        fa = self.write("p.fa", ">pep\nGFLG\n")
        out = self.tmp / "w.tsv"
        proc = self.run_script("protein_molecular_weight.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)["pep"]["Sequence"], "GFLG")

    def test_stdout_mode(self):
        fa = self.write("p.fa", ">pep\nGFLG\n")
        proc = self.run_script("protein_molecular_weight.py", "-i", fa)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Weight_Da", proc.stdout)
        self.assertIn("392.21", proc.stdout)

    def test_missing_input_clean_error(self):
        proc = self.run_script("protein_molecular_weight.py",
                               "-i", self.tmp / "no.fa")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)

    def test_empty_fasta_clean_error(self):
        fa = self.write("empty.fa", "")
        proc = self.run_script("protein_molecular_weight.py", "-i", fa)
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
