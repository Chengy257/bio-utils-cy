"""Functional tests for bin/fasta_to_alphafold_json.py."""

import json

from common import ScriptTestCase


class FastaToAlphafoldJsonTest(ScriptTestCase):
    def test_basic_conversion(self):
        fa = self.write("p.fa", ">prot1 some desc\nMKTAYIAKQRQ\n>prot2\nGFLG\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        data = json.loads(out.read_text())
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0]["name"], "prot1 some desc")
        self.assertEqual(data[0]["sequences"][0]["proteinChain"]["sequence"],
                         "MKTAYIAKQRQ")
        self.assertEqual(data[0]["sequences"][0]["proteinChain"]["count"], 1)
        self.assertEqual(data[0]["modelSeeds"], [])

    def test_lowercase_uppercased(self):
        fa = self.write("p.fa", ">prot\nmktayia\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        data = json.loads(out.read_text())
        self.assertEqual(data[0]["sequences"][0]["proteinChain"]["sequence"],
                         "MKTAYIA")

    def test_invalid_characters_skipped(self):
        fa = self.write("p.fa", ">good\nMKTAY\n>bad\nMK1*\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("bad", proc.stderr)
        data = json.loads(out.read_text())
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]["name"], "good")

    def test_empty_sequence_skipped(self):
        fa = self.write("p.fa", ">empty1\n>good\nMKT\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("empty sequence", proc.stderr)
        data = json.loads(out.read_text())
        self.assertEqual(len(data), 1)

    def test_all_invalid_exit_1(self):
        fa = self.write("p.fa", ">bad\n123\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 1)

    def test_chain_count(self):
        fa = self.write("p.fa", ">prot\nMKT\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out,
                               "--chain-count", "3")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        data = json.loads(out.read_text())
        self.assertEqual(data[0]["sequences"][0]["proteinChain"]["count"], 3)

    def test_invalid_chain_count_clean_error(self):
        fa = self.write("p.fa", ">prot\nMKT\n")
        out = self.tmp / "af.json"
        proc = self.run_script("fasta_to_alphafold_json.py", "-i", fa, "-o", out,
                               "--chain-count", "0")
        self.assertEqual(proc.returncode, 2)
        self.assertNotIn("Traceback", proc.stderr)

    def test_missing_input_clean_error(self):
        proc = self.run_script("fasta_to_alphafold_json.py",
                               "-i", self.tmp / "no.fa", "-o", self.tmp / "x.json")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
