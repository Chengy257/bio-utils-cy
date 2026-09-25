"""Functional tests for bin/kozak_similarity_score.py."""

from common import ScriptTestCase

# Per-position argmax bases of KOZAK_WEIGHTS (ATG rows are all-zero):
# score of this sequence must be exactly 1.0.
OPTIMAL = "CGCCGCCACCATGGCGGCGGAGG"
EXAMPLE = "GCCACCATGGCGATCGATCGATC"


class KozakSimilarityScoreTest(ScriptTestCase):
    def test_example_sequence_scores(self):
        proc = self.run_script("kozak_similarity_score.py", "-s", EXAMPLE)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Kozak similarity score:", proc.stdout)
        value = float(proc.stdout.strip().split(":")[1])
        self.assertGreaterEqual(value, 0.0)
        self.assertLessEqual(value, 1.0)

    def test_optimal_sequence_scores_one(self):
        proc = self.run_script("kozak_similarity_score.py", "-s", OPTIMAL)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        value = float(proc.stdout.strip().split(":")[1])
        self.assertAlmostEqual(value, 1.0, places=6)

    def test_wrong_length_clean_error(self):
        for seq in (OPTIMAL + "A", OPTIMAL[:-1]):
            proc = self.run_script("kozak_similarity_score.py", "-s", seq)
            self.assertEqual(proc.returncode, 1)
            self.assertNotIn("Traceback", proc.stderr)

    def test_non_atg_center_warns_but_scores(self):
        shifted = OPTIMAL[:10] + "ACG" + OPTIMAL[13:]
        proc = self.run_script("kozak_similarity_score.py", "-s", shifted)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("does not have ATG", proc.stderr)

    def test_rna_u_converted(self):
        seq = OPTIMAL[:10] + "AUG" + OPTIMAL[13:]
        proc = self.run_script("kozak_similarity_score.py", "-s", seq)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertNotIn("does not have ATG", proc.stderr)

    def test_batch_fasta(self):
        fa = self.write(
            "tis.fa",
            f">good\n{EXAMPLE}\n>bad_length\n{EXAMPLE}A\n>also_good\n{OPTIMAL}\n",
        )
        out = self.tmp / "scores.tsv"
        proc = self.run_script("kozak_similarity_score.py", "-i", fa, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Skipping >bad_length", proc.stderr)
        rows = out.read_text().strip().split("\n")
        self.assertEqual(rows[0], "seq_id\tKozak_score")
        self.assertEqual(len(rows), 3)  # header + 2 scored
        self.assertTrue(rows[1].startswith(">good\t"))
        self.assertEqual(rows[2].split("\t")[1], "1.000000")


if __name__ == "__main__":
    unittest.main()
