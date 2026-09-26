"""Functional tests for bin/aligned_fasta_similarity.py.

Hand-computed expectations over small synthetic alignments.
"""

import csv
import unittest

from common import ScriptTestCase

# 10 columns: M M -(-) M G(single) M M M M
# ref:    A C - A T - A G C A
# mouse:  A C - A A T A G C G
# -> double-gap col3 removed; col5 ref '-' (single gap), col6 mouse '-'
# comparable columns = 8, matches = 6 (cols 1,2,4,7,8,9) -> 0.75
ALN_MAIN = ">human.g1\nAC-AT-AGCA\n>mouse.g1\nAC-AATAGCG\n"

# ragged lengths (malformed alignment)
ALN_RAGGED = ">human.g1\nAC-AT-AGCA\n>rat.g1\nAC-AT-AG\n"


class AlignedFastaSimilarityTest(ScriptTestCase):

    def _dir_and_run(self, files, *extra):
        indir = self.tmp / "aln"
        indir.mkdir()
        for name, text in files.items():
            (indir / name).write_text(text)
        out = self.tmp / "sim.csv"
        proc = self.run_script("aligned_fasta_similarity.py", "-i", indir, "-r", "human",
                               "-o", out, *extra)
        return proc, out

    def test_hand_computed_metrics(self):
        proc, out = self._dir_and_run({"g.aligned.fa": ALN_MAIN})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["target"], "mouse.g1")
        self.assertEqual(int(r["matches"]), 6)
        self.assertEqual(int(r["gaps"]), 1)       # col6 single gap (col3 double-gap removed)
        self.assertEqual(int(r["total"]), 9)      # after double-gap removal
        self.assertAlmostEqual(float(r["similarity"]), 0.75)

    def test_all_gap_pair_similarity_none(self):
        # Every column has a gap on one side: comparable columns = 0,
        # similarity is empty in the CSV.
        aln = ">human.g1\nAC-A\n>fish.g1\n--A-\n"
        proc, out = self._dir_and_run({"g.aligned.fa": aln})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(rows[0]["similarity"], "")

    def test_mixed_gap_pair(self):
        # One comparable base column ('A'/'A'): similarity 1.0.
        aln = ">human.g1\nAC-A\n>fish.g1\n--AA\n"
        proc, out = self._dir_and_run({"g.aligned.fa": aln})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(float(rows[0]["similarity"]), 1.0)

    def test_ragged_alignment_skipped_with_warning(self):
        # Regression: ragged records used to be silently truncated and
        # compared anyway. (The second human.g1 here is a distinct
        # record object and is still compared — by design.)
        proc, out = self._dir_and_run({"g.aligned.fa": ALN_MAIN + ALN_RAGGED})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("length mismatch", proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual({r["target"] for r in rows}, {"mouse.g1", "human.g1"})
        self.assertNotIn("rat.g1", {r["target"] for r in rows})  # ragged -> skipped

    def test_reference_record_skipped_by_identity_not_id(self):
        # Two records share the reference's ID; only the chosen reference
        # record is skipped, the duplicate is still compared.
        aln = ">human.g1\nAC-AT-AGCA\n>human.g1\nAC-AT-AGCA\n>mouse.g1\nAC-AATAGCG\n"
        proc, out = self._dir_and_run({"g.aligned.fa": aln})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 2)
        self.assertEqual({r["target"] for r in rows}, {"human.g1", "mouse.g1"})
        # identical duplicate: similarity 1.0
        dup_row = next(r for r in rows if r["target"] == "human.g1")
        self.assertEqual(float(dup_row["similarity"]), 1.0)

    def test_missing_reference_no_results(self):
        proc, out = self._dir_and_run({"g.aligned.fa": ">mouse.g1\nACAT\n"})
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("not found", proc.stderr)
        self.assertFalse(out.exists())

    def test_threads_validation_and_no_files(self):
        proc, _ = self._dir_and_run({"g.aligned.fa": ALN_MAIN}, "-t", "0")
        self.assertEqual(proc.returncode, 2)

        indir = self.tmp / "empty"
        indir.mkdir()
        (indir / "x.txt").write_text("x")
        proc = self.run_script("aligned_fasta_similarity.py", "-i", indir,
                               "-r", "human", "-o", self.tmp / "o.csv")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
