"""Functional tests for bin/codon_frequency.py."""

from common import ScriptTestCase


class CodonFrequencyTest(ScriptTestCase):
    def read_table(self, out):
        rows = out.read_text().strip().split("\n")
        return rows[0], {r.split("\t")[0]: r.split("\t")[1:] for r in rows[1:]}

    def test_basic_counting_and_frequency(self):
        # ATG ATG GGG TTT -> 4 codons: ATG 2/4=50%, GGG 25%, TTT 25%
        fa = self.write("cds.fa", ">s1\nATGATGGGGTTT\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        header, table = self.read_table(out)
        self.assertEqual(header, "Codon\tCount\tFrequency(%)")
        self.assertEqual(table["ATG"], ["2", "50.0000"])
        self.assertEqual(table["GGG"], ["1", "25.0000"])
        self.assertEqual(table["TTT"], ["1", "25.0000"])

    def test_truncate(self):
        # 7 bases -> 2 codons with --truncate, clean error without
        fa = self.write("cds.fa", ">s1\nATGGGCT\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out, "--truncate")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertEqual(sum(int(v[0]) for v in table.values()), 2)

    def test_skip_invalid_denominator_excludes_skipped(self):
        # ATG NNN GGG -> with --skip-invalid: ATG 1/2=50%, GGG 50%
        fa = self.write("cds.fa", ">s1\nATGNNNGGG\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out, "--skip-invalid")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertEqual(table["ATG"], ["1", "50.0000"])
        self.assertEqual(table["GGG"], ["1", "50.0000"])
        self.assertNotIn("NNN", table)

    def test_invalid_codon_clean_error(self):
        fa = self.write("cds.fa", ">s1\nATGNNNGGG\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)

    def test_error_names_offending_sequence(self):
        fa = self.write("cds.fa", ">good one\nATGATGATG\n>bad seq\nATGAT\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertIn(">bad seq", proc.stderr)

    def test_lowercase_normalized(self):
        fa = self.write("cds.fa", ">s1\natgatggggttt\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertEqual(table["ATG"], ["2", "50.0000"])

    def test_t_and_u_counted_distinctly(self):
        fa = self.write("cds.fa", ">s1\nTTAUUA\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, table = self.read_table(out)
        self.assertIn("TTA", table)
        self.assertIn("UUA", table)

    def test_sort_freq(self):
        fa = self.write("cds.fa", ">s1\nATGATGGGGTTT\n")
        out = self.tmp / "freq.tsv"
        proc = self.run_script("codon_frequency.py", "-i", fa, "-o", out, "--sort", "freq")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = out.read_text().strip().split("\n")[1:]
        counts = [int(r.split("\t")[1]) for r in rows]
        self.assertEqual(counts, sorted(counts, reverse=True))


if __name__ == "__main__":
    unittest.main()
