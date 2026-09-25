"""Functional tests for bin/reverse_complement.py."""

from common import ScriptTestCase


class ReverseComplementTest(ScriptTestCase):
    def test_revcomp_string(self):
        proc = self.run_script("reverse_complement.py", "-s", "ATCGATCG", "-m", "revcomp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout.strip(), "CGATCGAT")

    def test_comp_and_rev_string(self):
        proc = self.run_script("reverse_complement.py", "-s", "ATCG", "-m", "comp")
        self.assertEqual(proc.stdout.strip(), "TAGC")
        proc = self.run_script("reverse_complement.py", "-s", "ATCG", "-m", "rev")
        self.assertEqual(proc.stdout.strip(), "GCTA")

    def test_fasta_roundtrip_keeps_description(self):
        fa = self.write("in.fa", ">gene1 some description here\nATCGATCG\n>gene2\nAATT\n")
        out = self.tmp / "out.fa"
        proc = self.run_script("reverse_complement.py", "-i", fa, "-m", "revcomp", "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(
            out.read_text(),
            ">gene1 some description here\nCGATCGAT\n>gene2\nAATT\n",
        )

    def test_rna_u_input_no_crash(self):
        # regression: U used to raise KeyError
        proc = self.run_script("reverse_complement.py", "-s", "AUCG", "-m", "revcomp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout.strip(), "CGAT")

    def test_gap_input_no_crash(self):
        # regression: - and . used to raise KeyError
        proc = self.run_script("reverse_complement.py", "-s", "AT-CG.", "-m", "revcomp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout.strip(), ".CG-AT")

    def test_invalid_character_clean_error(self):
        proc = self.run_script("reverse_complement.py", "-s", "ATZCG", "-m", "revcomp")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertIn("Z", proc.stderr)

    def test_missing_input_clean_error(self):
        proc = self.run_script("reverse_complement.py", "-i", self.tmp / "nope.fa")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)

    def test_lowercase_preserved(self):
        proc = self.run_script("reverse_complement.py", "-s", "atcg", "-m", "revcomp")
        self.assertEqual(proc.stdout.strip(), "cgat")


if __name__ == "__main__":
    unittest.main()
