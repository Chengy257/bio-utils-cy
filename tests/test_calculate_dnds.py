"""Functional tests for bin/calculate_dnds.py.

Expected counts hand-derived with the standard codon table:
- ATC (Ile) -> ATT (Ile)   : synonymous
- ATC (Ile) -> CTC (Leu)   : non-synonymous
- ATC (Ile) -> ATG (Met)   : non-synonymous
"""

import csv
import unittest

from common import ScriptTestCase


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class CalculateDndsUnitTest(unittest.TestCase):
    """Direct unit tests of the counting function."""

    @staticmethod
    def _calc(s1, s2):
        import importlib
        import sys
        from common import BIN
        sys.path.insert(0, str(BIN))
        dnds = importlib.import_module("calculate_dnds")
        return dnds.calculate_dn_ds(s1, s2)

    def test_synonymous_change(self):
        dn, ds, ratio = self._calc("ATC", "ATT")
        self.assertEqual((dn, ds), (0, 1))
        self.assertEqual(ratio, 0.0)

    def test_nonsynonymous_change(self):
        dn, ds, ratio = self._calc("ATC" + "ATC", "CTC" + "ATT")
        self.assertEqual((dn, ds), (1, 1))
        self.assertAlmostEqual(ratio, 1.0)

    def test_identical_sequences_na(self):
        dn, ds, ratio = self._calc("ATCGTT" * 4, "ATCGTT" * 4)
        self.assertEqual((dn, ds, ratio), (0, 0, None))

    def test_all_nonsynonymous_inf(self):
        # ATC->CTC, ATC->CTA? (CTA=Leu too) — keep it simple: Ile->Leu both.
        dn, ds, ratio = self._calc("ATCATC", "CTCCTC")
        self.assertEqual((dn, ds), (2, 0))
        self.assertEqual(ratio, float("inf"))

    def test_gaps_and_partial_tail_ignored(self):
        # Gap codons skipped; trailing partial codon (2 nt) ignored.
        dn, ds, ratio = self._calc("ATC---ATTAC", "ATTCGGATTAC")
        self.assertEqual((dn, ds), (0, 1))

    def test_unequal_lengths_raise(self):
        with self.assertRaises(ValueError):
            self._calc("ATCG", "ATC")


class CalculateDndsTest(ScriptTestCase):

    def _dir_and_run(self, files):
        indir = self.tmp / "alignments"
        indir.mkdir()
        for name, text in files.items():
            (indir / name).write_text(text)
        out = self.tmp / "dnds.tsv"
        proc = self.run_script("calculate_dnds.py", "-i", indir, "-o", out)
        return proc, out

    def test_pair_counts_and_na_ratio(self):
        # sp1 vs sp2: one syn + one nonsyn -> dN=1 dS=1 ratio=1.0
        # sp3 identical to sp1 -> NA
        fasta = (
            ">sp1\nATCATCGGT\n"
            ">sp2\nCTCATTGGT\n"
            ">sp3\nATCATCGGT\n"
        )
        proc, out = self._dir_and_run({"gene1.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        by_pair = {(r["Species1"], r["Species2"]): r for r in rows}
        r12 = by_pair[("sp1", "sp2")]
        self.assertEqual((r12["dN"], r12["dS"], r12["dN_dS"]), ("1", "1", "1.0000"))
        r13 = by_pair[("sp1", "sp3")]
        self.assertEqual((r13["dN"], r13["dS"], r13["dN_dS"]), ("0", "0", "NA"))
        # Header check
        with open(out) as fh:
            self.assertEqual(fh.readline().strip(), "Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS")

    def test_lowercase_input_normalized(self):
        fasta = ">sp1\natcatcggt\n>sp2\nCTCATTGGT\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual((rows[0]["dN"], rows[0]["dS"]), ("1", "1"))

    def test_protein_length_mismatch_skips_pair_visibly(self):
        # sp2 carries an internal stop: TAA at codon 2 truncates its
        # protein to 1 aa vs 3 aa -> pair skipped, reported at INFO.
        fasta = (
            ">sp1\nATCATCATC\n"
            ">sp2\nATCTAATTC\n"
        )
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("length mismatch", proc.stderr)
        with open(out) as fh:
            body = fh.readlines()[1:]
        self.assertEqual(body, [])

    def test_trailing_partial_codon_ignored(self):
        fasta = ">sp1\nATCATCCG\n>sp2\nATCATTAA\n"  # last codon CG/AA invalid -> skipped
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        # 8 nt -> floor codons: ATC ATC / ATC ATT compared, 2-nt tail dropped
        self.assertEqual((rows[0]["dN"], rows[0]["dS"]), ("0", "1"))

    def test_duplicate_ids_warn(self):
        fasta = ">sp1\nATCATCGGT\n>sp1\nCTCATTGGT\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Duplicate sequence ID", proc.stderr)

    def test_batch_multiple_files_and_suffix(self):
        f1 = ">a\nATCATCGGT\n>b\nCTCATTGGT\n"
        f2 = ">a\nATCGTTTTT\n>b\nATCGTTTTA\n"  # syn change TTT->TTA (Phe)
        proc, out = self._dir_and_run({"g1.dna.fa": f1, "g2.dna.fa": f2, "notes.txt": "skip me"})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual({r["Gene"] for r in rows}, {"g1", "g2"})
        self.assertEqual(len(rows), 2)

    def test_no_matching_files_exits_1(self):
        proc, _ = self._dir_and_run({"readme.txt": "hello"})
        self.assertEqual(proc.returncode, 1)

    def test_custom_suffix(self):
        f1 = ">a\nATCATCGGT\n>b\nCTCATTGGT\n"
        proc, out = self._dir_and_run({"g1.fa": f1})
        self.assertEqual(proc.returncode, 1)  # default suffix misses
        indir = self.tmp / "alignments"
        out = self.tmp / "dnds.tsv"
        proc = self.run_script("calculate_dnds.py", "-i", indir, "-o", out, "--suffix", ".fa")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(read_tsv(out)), 1)


if __name__ == "__main__":
    unittest.main()
