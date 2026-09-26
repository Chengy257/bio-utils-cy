"""Functional tests for bin/blast_align_analysis.py.

blastn/blastp are on PATH, so these run real BLAST comparisons on
~120 nt synthetic CDS pairs.
"""

import csv
import unittest

from common import ScriptTestCase

DNA_HUMAN = (
    "ATGGCCTCTGGCCCCCTGGTGCCCTGCCTGGCTTCCTGCCATTCTGGTGCTGGTGGCCCTGG"
    "GCTCCTGCCTGCTGGCTTCCTGCCATTCTGGTGCTGGTGGCCCTGGGCTCCTGCCTGCTGGC"
)
# Same 120 nt with a handful of quiet substitutions (~90% identity).
DNA_MOUSE = (
    "ATGGCCTCTGGCCCCCTGGTGCCCTGCCCAGCTTCCTGCCATTCTGGTGCTGGTGGCCCTGG"
    "GCTCCTGCGTGCTGGCTTCCTGCCATTCTGGTGCTGGTGGCCCTGGGCTCCTGCCTGCTGA "
).replace(" ", "C")
# No meaningful similarity to the others (different composition/frame).
DNA_ZEBRA = "ACACACACAC" * 5 + "GTGTGTGTGT" * 5 + "ACACACACAC" * 2


class BlastAlignAnalysisTest(ScriptTestCase):

    def _dir_and_run(self, fasta_text, *extra):
        indir = self.tmp / "aln"
        indir.mkdir()
        (indir / "gene1.dna.fa").write_text(fasta_text)
        out = self.tmp / "res.csv"
        proc = self.run_script("blast_align_analysis.py", "-i", indir, "-r", "human",
                               "-o", out, "-t", "1", *extra)
        return proc, out

    def test_close_pair_high_identity(self):
        fasta = f">human.g1\n{DNA_HUMAN}\n>mouse.g1\n{DNA_MOUSE}\n>zebra.g1\n{DNA_ZEBRA}\n"
        proc, out = self._dir_and_run(fasta)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = {r["target"]: r for r in csv.DictReader(fh)}
        self.assertEqual(set(rows), {"mouse.g1", "zebra.g1"})
        m = rows["mouse.g1"]
        self.assertIsNotNone(m["nucleotide_identity"])
        self.assertGreater(float(m["nucleotide_identity"]), 0.85)
        self.assertLess(float(m["nucleotide_evalue"]), 1e-10)
        self.assertIsNotNone(m["protein_identity"])
        # Unrelated pair: no BLAST alignment -> empty cells.
        z = rows["zebra.g1"]
        self.assertTrue(z["nucleotide_evalue"] in ("", "None") or float(z["nucleotide_evalue"]) > 1e-3)

    def test_custom_task_blastn(self):
        # The standard blastn task (11-mer seeds) cannot seed a 124 nt,
        # 4-mismatch pair at all; give it an exact 60 nt stretch so a
        # seed exists. This pins the --task plumbing end to end.
        shared = DNA_HUMAN[:60]
        divergent = DNA_ZEBRA[:60]
        fasta = f">human.g1\n{shared}{DNA_HUMAN[60:]}\n>mouse.g1\n{shared}{divergent}\n"
        proc, out = self._dir_and_run(fasta, "--task", "blastn")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 1)
        self.assertGreaterEqual(float(rows[0]["nucleotide_identity"]), 0.5)

    def test_invalid_task_rejected(self):
        fasta = f">human.g1\n{DNA_HUMAN}\n>mouse.g1\n{DNA_MOUSE}\n"
        proc, _ = self._dir_and_run(fasta, "--task", "tblastx")
        self.assertEqual(proc.returncode, 2)

    def test_missing_reference_species(self):
        fasta = f">mouse.g1\n{DNA_MOUSE}\n"
        proc, out = self._dir_and_run(fasta)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("not found", proc.stderr)
        self.assertFalse(out.exists())  # no results -> no CSV

    def test_internal_stop_codon_does_not_reach_blastp(self):
        # TAA inside human.g1 at codon 2: to_stop translation truncates
        # instead of feeding '*' to blastp (which used to break).
        stop_mid = DNA_HUMAN[:6] + "TAA" + DNA_HUMAN[9:]
        fasta = f">human.g1\n{stop_mid}\n>mouse.g1\n{DNA_MOUSE}\n"
        proc, out = self._dir_and_run(fasta)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue(out.exists())

    def test_threads_validation_and_no_files(self):
        fasta = f">human.g1\n{DNA_HUMAN}\n>mouse.g1\n{DNA_MOUSE}\n"
        proc, _ = self._dir_and_run(fasta, "-t", "0")
        self.assertEqual(proc.returncode, 2)

        indir = self.tmp / "empty_aln"
        indir.mkdir()
        (indir / "notes.txt").write_text("x")
        proc = self.run_script("blast_align_analysis.py", "-i", indir, "-r", "human",
                               "-o", self.tmp / "o.csv")
        self.assertEqual(proc.returncode, 1)

    def test_gap_only_sequence_fails_validation(self):
        fasta = ">human.g1\nATGGCCTCTGGCCCCCTGGTGCCCTGCCTGGC\n>bad.g1\n----\n"
        proc, _ = self._dir_and_run(fasta)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Invalid sequence", proc.stderr)


if __name__ == "__main__":
    unittest.main()
