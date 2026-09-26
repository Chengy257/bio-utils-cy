"""Functional tests for bin/sequence_similarity_network.py.

python-louvain (community) is not installed in this environment, so the
Louvain step cannot run here; everything else (k-mer profile, binary
Jaccard matrix, network build, SIF, PDF with a stubbed partition, and
the clean missing-dependency abort) is tested for real.
"""

import importlib
import sys
import unittest
from pathlib import Path

from common import BIN, ScriptTestCase

sys.path.insert(0, str(BIN))
ssn = importlib.import_module("sequence_similarity_network")


def fasta(paths_seqs):
    return "".join(f">{n}\n{s}\n" for n, s in paths_seqs)


class SequenceSimilarityNetworkUnitTest(unittest.TestCase):

    def test_get_kmers(self):
        self.assertEqual(ssn.get_kmers("ACGTA", 3), ["ACG", "CGT", "GTA"])
        self.assertEqual(ssn.get_kmers("AC", 3), [])

    def test_identical_sequences_full_similarity(self):
        seqs = ["ACGTACGTA", "ACGTACGTA"]
        m = ssn.compute_similarity_matrix(seqs, k=3, threshold=0.0)
        self.assertEqual(round(m[0, 1], 6), 1.0)

    def test_binary_jaccard_semantics(self):
        # Binary presence/absence: two sequences with identical k-mer
        # SETS score 1.0 even when k-mer counts differ wildly.
        # "AAATAAATAAAT" is period-3, so doubling introduces no new
        # 3-mers at the junction.
        s1 = "AAATAAATAAAT" * 2
        s2 = "AAATAAAT"
        m = ssn.compute_similarity_matrix([s1, s2], k=3, threshold=0.0)
        self.assertEqual(round(m[0, 1], 6), 1.0)

    def test_disjoint_sequences_below_threshold_zeroed(self):
        s1 = "ACACACACAC"
        s2 = "TGTGTGTGTG"
        m = ssn.compute_similarity_matrix([s1, s2], k=3, threshold=0.5)
        self.assertEqual(m[0, 1], 0)

    def test_fewer_than_two_sequences_exits(self):
        with self.assertRaises(SystemExit) as ctx:
            ssn.compute_similarity_matrix(["ACGT"], k=3)
        self.assertEqual(ctx.exception.code, 1)


class SequenceSimilarityNetworkTest(ScriptTestCase):

    def _fasta(self, name, pairs):
        return self.write(name, fasta(pairs))

    def test_help_and_version_work_without_dependencies(self):
        proc = self.run_script("sequence_similarity_network.py", "--help")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        proc = self.run_script("sequence_similarity_network.py", "--version")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("1.1.0", proc.stdout)

    def test_missing_python_louvain_clean_exit_1(self):
        # community is not installed: run up to the Louvain step and
        # assert the clean abort (no traceback, exit 1).
        inp = self._fasta("in.fa", [
            ("a1", "ACGTACGTAC" * 3),
            ("a2", "ACGTACGTAC" * 2 + "ACGTAC"),  # shares k-mers with a1
            ("b1", "TTTTGGGGCC" * 3),
            ("b2", "TTTTGGGGCC" * 2 + "TTTGG"),   # shares k-mers with b1
        ])
        proc = self.run_script("sequence_similarity_network.py", "-i", inp,
                               "-o", self.tmp / "out", "-k", "3", "-t", "0.5")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("python-louvain", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_threshold_validation(self):
        inp = self._fasta("in.fa", [("a", "ACGTACGTAC"), ("b", "ACGTACGTAC")])
        proc = self.run_script("sequence_similarity_network.py", "-i", inp,
                               "-o", self.tmp / "out", "-t", "1.5")
        self.assertEqual(proc.returncode, 2)  # argparse error
        self.assertIn("threshold", proc.stderr)

    def test_kmer_validation(self):
        inp = self._fasta("in.fa", [("a", "ACGT"), ("b", "ACGT")])
        proc = self.run_script("sequence_similarity_network.py", "-i", inp,
                               "-o", self.tmp / "out", "-k", "0")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("kmer-size", proc.stderr)

    def test_missing_or_empty_input(self):
        proc = self.run_script("sequence_similarity_network.py", "-i",
                               self.tmp / "nope.fa", "-o", self.tmp / "out")
        self.assertEqual(proc.returncode, 1)
        empty = self.write("empty.fa", "")
        proc = self.run_script("sequence_similarity_network.py", "-i", empty,
                               "-o", self.tmp / "out")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No sequences", proc.stderr)

    def test_pipeline_with_stubbed_partition(self):
        # Stub Louvain with a fixed 2-community partition; the rest of
        # the pipeline (matrix, network, PDF, SIF) runs for real.
        pairs = [
            ("a1", "ACGTACGTAC" * 4),
            ("a2", "ACGTACGTAC" * 3 + "ACGTAC"),
            ("a3", "ACGTACGTAC" * 2 + "ACGTACGTAC"),
            ("b1", "TTTTGGGGCC" * 4),
            ("b2", "TTTTGGGGCC" * 3 + "TTTGG"),
        ]
        inp = self._fasta("in.fa", pairs)
        prefix = self.tmp / "net"

        real_detect = ssn.detect_communities
        ssn.detect_communities = lambda G: {n: (0 if n < 3 else 1) for n in G.nodes()}
        # Route main() to the temp dir for outputs.
        old_argv = sys.argv
        try:
            sys.argv = [
                "sequence_similarity_network.py", "-i", str(inp),
                "-o", str(prefix), "-k", "4", "-t", "0.5",
            ]
            ssn.main()
        finally:
            sys.argv = old_argv
            ssn.detect_communities = real_detect

        sif = Path(f"{prefix}_network.pdf".replace("_network.pdf", "_network.sif"))
        pdf = Path(f"{prefix}_network.pdf")
        self.assertTrue(pdf.exists() and pdf.stat().st_size > 1000)
        lines = sif.read_text().splitlines()
        self.assertTrue(lines)
        for ln in lines:
            a, _, b = ln.split("\t")
            self.assertIn({a, b}, [{"a1", "a2"}, {"a1", "a3"}, {"a2", "a3"},
                                   {"b1", "b2"}])
        # a-cluster must not connect to b-cluster.
        for ln in lines:
            self.assertNotIn("b1\tsimilarity\ta", ln)


if __name__ == "__main__":
    unittest.main()
