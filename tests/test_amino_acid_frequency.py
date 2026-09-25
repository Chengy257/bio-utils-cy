"""Functional tests for bin/amino_acid_frequency.py."""

from common import ScriptTestCase


class AminoAcidFrequencyTest(ScriptTestCase):
    def read_table(self, out):
        rows = out.read_text().strip().split("\n")
        return rows[0], {r.split("\t")[0]: r.split("\t")[1:] for r in rows[1:]}

    def test_basic_counts(self):
        # PEPTIDEK: P1 E2 T1 I1 D1 K1 -> 8 residues
        fa = self.write("p.fa", ">t1\nPEPTIDEK\n")
        out = self.tmp / "aa.tsv"
        proc = self.run_script("amino_acid_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        header, table = self.read_table(out)
        self.assertEqual(header, "Amino_acid\tCount\tFrequency(%)")
        self.assertEqual(table["E"], ["2", "25.0000"])
        self.assertEqual(table["K"], ["1", "12.5000"])

    def test_standard_only_denominator_excludes_nonstandard(self):
        # PEPTIDEKX* : 8 standard + X + * -> standard-only: E 2/8=25%, sum 100
        fa = self.write("p.fa", ">t1\nPEPTIDEKX*\n")
        out = self.tmp / "aa.tsv"
        proc = self.run_script("amino_acid_frequency.py", "-i", fa, "-o", out, "--standard-only")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertNotIn("X", table)
        self.assertNotIn("*", table)
        total_freq = sum(float(v[1]) for v in table.values())
        self.assertAlmostEqual(total_freq, 100.0, places=3)
        self.assertEqual(table["E"], ["2", "25.0000"])

    def test_default_counts_every_character(self):
        fa = self.write("p.fa", ">t1\nPEPX\n")
        out = self.tmp / "aa.tsv"
        proc = self.run_script("amino_acid_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertIn("X", table)
        self.assertEqual(table["X"], ["1", "25.0000"])

    def test_lowercase_normalized(self):
        fa = self.write("p.fa", ">t1\npeptk\n")
        out = self.tmp / "aa.tsv"
        proc = self.run_script("amino_acid_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertEqual(table["P"], ["2", "40.0000"])
        self.assertEqual(table["K"], ["1", "20.0000"])
        self.assertNotIn("p", table)

    def test_empty_file_clean_error(self):
        fa = self.write("p.fa", ">header_only\n")
        out = self.tmp / "aa.tsv"
        proc = self.run_script("amino_acid_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)

    def test_sort_freq(self):
        fa = self.write("p.fa", ">t1\nPEPTIDEK\n")
        out = self.tmp / "aa.tsv"
        proc = self.run_script("amino_acid_frequency.py", "-i", fa, "-o", out, "--sort", "freq")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = out.read_text().strip().split("\n")[1:]
        counts = [int(r.split("\t")[1]) for r in rows]
        self.assertEqual(counts, sorted(counts, reverse=True))


if __name__ == "__main__":
    unittest.main()
