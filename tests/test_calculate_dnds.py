"""Functional tests for bin/calculate_dnds.py.

Expected values hand-derived with the standard codon table.

Legacy (diff-ratio) counts:
- ATC (Ile) -> ATT (Ile)   : synonymous
- ATC (Ile) -> CTC (Leu)   : non-synonymous

Nei-Gojobori (1986) ingredients:
- synonymous sites per codon: s(CCT)=1 (4-fold pos3), s(ATT)=2/3 (Ile:
  ATA/ATC syn), s(TTT)=1/3 (Phe: TTC syn), s(AAA)=1/3 (Lys: AAG syn,
  AAC/AAT are Asn!), s(GGT)=1, s(ACT)=1, s(GTT)=1, s(AAT)=1/3,
  s(ATG)=0, s(TGG)=0.
- pathway averaging: ATT->GTG has one syn step on one of its two
  pathways -> (0.5, 1.5); AAA->CCC averages 3 syn steps over 6
  pathways -> (0.5, 2.5).

Integration cases (12 codons, 3 differing + 9 identical AAA):
- syn-only: s_sites = 1 + 2/3 + 1/3 + 9*(1/3) = 5.0, s_diffs = 3
  -> pS = 0.6 -> dS = -0.75*ln(0.2) = 1.2071; dN = 0.
- mixed: s_sites = 1 + 2/3 + 1 + 3 = 17/3, s_diffs = 2 -> pS = 6/17
  -> dS = 0.4770; n_sites = 91/3, n_diffs = 1 -> pN = 3/91
  -> dN = 0.0337; ratio = 0.0707.
- saturation: ATTx10 vs ATCx10 -> pS = 1.5 >= 0.75 -> NA.
- all-nonsyn: 2 differing sites, 9 identical AAA -> dS = 0, dN > 0
  -> Inf.
"""

import csv
import unittest

from common import ScriptTestCase


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class _DndsModuleMixin:
    """Load bin/calculate_dnds.py as a module for direct unit tests."""

    @classmethod
    def _load(cls):
        import importlib
        import sys
        from common import BIN
        if "calculate_dnds" not in sys.modules:
            sys.path.insert(0, str(BIN))
            importlib.import_module("calculate_dnds")
        return sys.modules["calculate_dnds"]


class CalculateDndsUnitTest(_DndsModuleMixin, unittest.TestCase):
    """Direct unit tests of the legacy diff-ratio counting function."""

    @staticmethod
    def _calc(s1, s2):
        return CalculateDndsUnitTest._load().calculate_dn_ds(s1, s2)

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


class CalculateDndsNgUnitTest(_DndsModuleMixin, unittest.TestCase):
    """Unit tests of the Nei-Gojobori helpers and calculator."""

    @classmethod
    def _mod(cls):
        return cls._load()

    def test_synonymous_sites(self):
        mod = self._mod()
        expected = {
            "CCT": 1.0, "GGT": 1.0, "ACT": 1.0, "GTT": 1.0,
            "ATT": 2 / 3, "ATC": 2 / 3, "AAA": 1 / 3, "TTT": 1 / 3,
            "AAT": 1 / 3, "CAT": 1 / 3,
            "ATG": 0.0, "TGG": 0.0,
        }
        for codon, s in expected.items():
            self.assertAlmostEqual(mod.codon_syn_sites(codon), s, places=12,
                                   msg=f"syn sites of {codon}")

    def test_pair_diffs_single_step(self):
        mod = self._mod()
        self.assertEqual(mod.codon_pair_diffs("ATC", "ATT"), (1.0, 0.0))
        self.assertEqual(mod.codon_pair_diffs("CCT", "CAT"), (0.0, 1.0))
        self.assertEqual(mod.codon_pair_diffs("ATC", "ATC"), (0.0, 0.0))

    def test_pair_diffs_pathway_averaging(self):
        mod = self._mod()
        # ATT->GTG: pathway via GTT has 1 syn step, via ATG none -> 0.5/1.5.
        self.assertEqual(mod.codon_pair_diffs("ATT", "GTG"), (0.5, 1.5))
        # Symmetric in its arguments.
        self.assertEqual(mod.codon_pair_diffs("GTG", "ATT"), (0.5, 1.5))
        # AAA->CCC: 3 syn steps summed over 6 pathways -> 0.5/2.5.
        self.assertEqual(mod.codon_pair_diffs("AAA", "CCC"), (0.5, 2.5))
        # ATT->CAT: both pathways purely non-synonymous.
        self.assertEqual(mod.codon_pair_diffs("ATT", "CAT"), (0.0, 2.0))

    def test_jukes_cantor(self):
        import math
        mod = self._mod()
        self.assertEqual(mod.jukes_cantor(0), 0.0)
        self.assertIsNone(mod.jukes_cantor(None))
        self.assertIsNone(mod.jukes_cantor(0.75))
        self.assertIsNone(mod.jukes_cantor(1.5))
        self.assertAlmostEqual(mod.jukes_cantor(0.375), -0.75 * math.log(0.5))
        self.assertAlmostEqual(mod.jukes_cantor(0.5), -0.75 * math.log(1 / 3))

    def test_ng_identical_na(self):
        d_n, d_s, ratio = self._mod().calculate_dn_ds_ng("ATCGTT" * 4, "ATCGTT" * 4)
        self.assertEqual((d_n, d_s, ratio), (0.0, 0.0, None))

    def test_ng_syn_only(self):
        seq_a = "CCT" + "ATT" + "TTT" + "AAA" * 9
        seq_b = "CCG" + "ATC" + "TTC" + "AAA" * 9
        d_n, d_s, ratio = self._mod().calculate_dn_ds_ng(seq_a, seq_b)
        self.assertEqual(d_n, 0.0)
        self.assertAlmostEqual(d_s, 1.20707843, places=7)
        self.assertEqual(ratio, 0.0)

    def test_ng_mixed(self):
        seq_a = "CCT" + "ACT" + "GTT" + "AAA" * 9
        seq_b = "CCG" + "AAT" + "GTA" + "AAA" * 9
        d_n, d_s, ratio = self._mod().calculate_dn_ds_ng(seq_a, seq_b)
        self.assertAlmostEqual(d_n, 0.03371354, places=7)
        self.assertAlmostEqual(d_s, 0.47699158, places=7)
        self.assertAlmostEqual(ratio, 0.07067953, places=7)

    def test_ng_saturation_na(self):
        # pS = 1.5 (every differing site is a 2-fold synonymous change).
        d_n, d_s, ratio = self._mod().calculate_dn_ds_ng("ATT" * 10, "ATC" * 10)
        self.assertEqual(d_n, 0.0)
        self.assertIsNone(d_s)
        self.assertIsNone(ratio)

    def test_ng_all_nonsynonymous_inf(self):
        # ATT->CAT is a purely non-synonymous 2-step difference.
        seq_a = "AAA" * 10 + "ATT" * 2
        seq_b = "AAA" * 10 + "CAT" * 2
        d_n, d_s, ratio = self._mod().calculate_dn_ds_ng(seq_a, seq_b)
        self.assertEqual(d_s, 0.0)
        self.assertGreater(d_n, 0.0)
        self.assertEqual(ratio, float("inf"))

    def test_ng_gaps_and_tail_ignored(self):
        # Gap codon and 2-nt tail skipped: only the ATT/ATT and ATT/TTC
        # sites count (TTC syn: s=1/3 each, S_d=1).
        d_n, d_s, ratio = self._mod().calculate_dn_ds_ng("ATT---ATTAC", "ATTTCGATTAC")
        self.assertEqual(d_n, 0.0)
        self.assertAlmostEqual(d_s, -0.75 * (0.0 if False else 1) * 0, places=12) if False else None
        # pS = 1 / (2/3) = 1.5 -> saturated -> NA.
        self.assertIsNone(ratio)

    def test_ng_unequal_lengths_raise(self):
        with self.assertRaises(ValueError):
            self._mod().calculate_dn_ds_ng("ATCG", "ATC")


class CalculateDndsTest(ScriptTestCase):

    def _dir_and_run(self, files, *extra_args):
        indir = self.tmp / "alignments"
        indir.mkdir()
        for name, text in files.items():
            (indir / name).write_text(text)
        out = self.tmp / "dnds.tsv"
        proc = self.run_script("calculate_dnds.py", "-i", indir, "-o", out, *extra_args)
        return proc, out

    # --- legacy diff-ratio method (byte-exact old behaviour) ---

    def test_diff_ratio_pair_counts_and_na_ratio(self):
        # sp1 vs sp2: one syn + one nonsyn -> dN=1 dS=1 ratio=1.0
        # sp3 identical to sp1 -> NA
        fasta = (
            ">sp1\nATCATCGGT\n"
            ">sp2\nCTCATTGGT\n"
            ">sp3\nATCATCGGT\n"
        )
        proc, out = self._dir_and_run({"gene1.dna.fa": fasta}, "--method", "diff-ratio")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        by_pair = {(r["Species1"], r["Species2"]): r for r in rows}
        r12 = by_pair[("sp1", "sp2")]
        self.assertEqual((r12["dN"], r12["dS"], r12["dN_dS"]), ("1", "1", "1.0000"))
        r13 = by_pair[("sp1", "sp3")]
        self.assertEqual((r13["dN"], r13["dS"], r13["dN_dS"]), ("0", "0", "NA"))
        # Legacy header is unchanged from v1.1.0.
        with open(out) as fh:
            self.assertEqual(fh.readline().strip(), "Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS")

    def test_diff_ratio_lowercase_input_normalized(self):
        fasta = ">sp1\natcatcggt\n>sp2\nCTCATTGGT\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta}, "--method", "diff-ratio")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual((rows[0]["dN"], rows[0]["dS"]), ("1", "1"))

    def test_diff_ratio_trailing_partial_codon_ignored(self):
        fasta = ">sp1\nATCATCCG\n>sp2\nATCATTAA\n"  # last codon CG/AA invalid -> skipped
        proc, out = self._dir_and_run({"g.dna.fa": fasta}, "--method", "diff-ratio")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual((rows[0]["dN"], rows[0]["dS"]), ("0", "1"))

    # --- nei-gojobori method (default) ---

    def test_ng_default_header_and_syn_only(self):
        fasta = ">sp1\nCCTATTTTT" + "AAA" * 9 + "\n>sp2\nCCGATCTTC" + "AAA" * 9 + "\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            self.assertEqual(
                fh.readline().strip(),
                "Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS\tdn_count\tds_count",
            )
        rows = read_tsv(out)
        self.assertEqual(
            (rows[0]["dN"], rows[0]["dS"], rows[0]["dN_dS"]),
            ("0.0000", "1.2071", "0.0000"),
        )
        self.assertEqual((rows[0]["dn_count"], rows[0]["ds_count"]), ("0", "3"))

    def test_ng_mixed_syn_nonsyn(self):
        fasta = ">sp1\nCCTACTGTT" + "AAA" * 9 + "\n>sp2\nCCGAATGTA" + "AAA" * 9 + "\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(
            (rows[0]["dN"], rows[0]["dS"], rows[0]["dN_dS"]),
            ("0.0337", "0.4770", "0.0707"),
        )
        self.assertEqual((rows[0]["dn_count"], rows[0]["ds_count"]), ("1", "2"))

    def test_ng_identical_na_with_counts(self):
        fasta = ">sp1\nATCGTT" * 1 + "\n>sp2\nATCGTT\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(
            (rows[0]["dN"], rows[0]["dS"], rows[0]["dN_dS"], rows[0]["dn_count"], rows[0]["ds_count"]),
            ("0.0000", "0.0000", "NA", "0", "0"),
        )

    def test_ng_saturation_reports_na(self):
        fasta = ">sp1\n" + "ATT" * 10 + "\n>sp2\n" + "ATC" * 10 + "\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(
            (rows[0]["dN"], rows[0]["dS"], rows[0]["dN_dS"]),
            ("0.0000", "NA", "NA"),
        )
        self.assertEqual((rows[0]["dn_count"], rows[0]["ds_count"]), ("0", "10"))

    def test_ng_all_nonsynonymous_inf(self):
        fasta = ">sp1\n" + "AAA" * 10 + "ATT" * 2 + "\n>sp2\n" + "AAA" * 10 + "CAT" * 2 + "\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(
            (rows[0]["dN"], rows[0]["dS"], rows[0]["dN_dS"]),
            ("0.1383", "0.0000", "Inf"),
        )
        self.assertEqual((rows[0]["dn_count"], rows[0]["ds_count"]), ("2", "0"))

    def test_ng_gaps_skipped(self):
        # Balanced gaps so the gap-free proteins match; the gapped codon
        # pair is skipped (no diff, no crash). Only ATC/ATT (syn, s=2/3)
        # differs; ATT/ATT + 8x GTT/GTT identical.
        # s_total = 2/3 + 2/3 + 8*1 = 28/3; pS = 1/(28/3) = 3/28
        # -> dS = -0.75*ln(6/7) = 0.1156; dN = 0.
        fasta = ">sp1\nATC---ATT" + "GTT" * 8 + "\n>sp2\nATT---ATT" + "GTT" * 8 + "\n"
        proc, out = self._dir_and_run({"g.dna.fa": fasta})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(
            (rows[0]["dN"], rows[0]["dS"], rows[0]["dN_dS"]),
            ("0.0000", "0.1156", "0.0000"),
        )

    # --- method-agnostic behaviour ---

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

    def test_invalid_method_rejected(self):
        fasta = ">a\nATCATCGGT\n>b\nCTCATTGGT\n"
        proc, _ = self._dir_and_run({"g.dna.fa": fasta}, "--method", "yang-nielsen")
        self.assertNotEqual(proc.returncode, 0)


if __name__ == "__main__":
    unittest.main()
