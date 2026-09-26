"""Functional tests for bin/blast_sequence_network.py.

markov-clustering is not installed in this environment; the MCL step is
exercised with a stub module injected into sys.modules, while BLAST
parsing, the real blastp all-vs-all run, and the visualization are
tested for real.
"""

import importlib
import sys
import types
import unittest
from pathlib import Path

from common import BIN, ScriptTestCase

sys.path.insert(0, str(BIN))
bsn = importlib.import_module("blast_sequence_network")


def protein_fa(pairs):
    return "".join(f">{n}\n{s}\n" for n, s in pairs)


# Two families of clearly-related proteins.
FAM_A = ("a1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFVVQLFE")
FAM_A2 = ("a2", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFVVQLFE")
FAM_A3 = ("a3", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFVVQLFD")
FAM_B = ("b1", "GSDQWNKELGAQALEHHSQIGLLYTPGQQLAFHASNNQTLGTL FAGGQYQELMAAQVQSPQDFSLILQ")
FAM_B = ("b1", FAM_B[1].replace(" ", ""))
FAM_B2 = ("b2", FAM_B[1][:-4] + "MAAQ")


class BlastSequenceNetworkUnitTest(unittest.TestCase):

    def test_parse_skips_self_hits_and_applies_threshold(self):
        blast_txt = (
            "a1\ta1\t100.000\t80\n"   # self-hit: must be ignored
            "a1\ta2\t95.000\t80\n"    # kept (0.95 >= 0.5)
            "a1\tb1\t30.000\t80\n"    # dropped (< 0.5)
            "a2\ta1\t90.000\t80\n"    # symmetric fill
            "b1\tb2\t99.000\t80\n"
        )
        p = Path(self.tmpdir()) / "b.txt"
        p.write_text(blast_txt)
        matrix, ids = bsn.parse_blast_results(str(p), ["a1", "a2", "b1", "b2"], 0.5)
        self.assertEqual(ids, ["a1", "a2", "b1", "b2"])
        self.assertEqual(matrix[0, 0], 0)        # diagonal zeroed
        # Symmetric fill; the later (a2,a1) line overwrites (last wins).
        self.assertEqual(round(matrix[0, 1], 2), 0.90)
        self.assertEqual(round(matrix[1, 0], 2), 0.90)
        self.assertEqual(matrix[0, 2], 0)
        self.assertEqual(round(matrix[2, 3], 2), 0.99)

    def test_threshold_scale_one_half_means_fifty_percent(self):
        blast_txt = "a1\ta2\t50.000\t80\n"
        p = Path(self.tmpdir()) / "b.txt"
        p.write_text(blast_txt)
        # 0.5 keeps a 50% hit; 0.6 does not.
        m, _ = bsn.parse_blast_results(str(p), ["a1", "a2"], 0.5)
        self.assertEqual(round(m[0, 1], 2), 0.50)
        m, _ = bsn.parse_blast_results(str(p), ["a1", "a2"], 0.6)
        self.assertEqual(m[0, 1], 0)

    def tmpdir(self):
        import tempfile
        return tempfile.mkdtemp()


class BlastSequenceNetworkTest(ScriptTestCase):

    def _fasta(self, name="in.fa"):
        return self.write(name, protein_fa([FAM_A, FAM_A2, FAM_A3, FAM_B, FAM_B2]))

    def test_help_version_identity_threshold_validation(self):
        proc = self.run_script("blast_sequence_network.py", "--help")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        proc = self.run_script("blast_sequence_network.py", "--version")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("1.1.0", proc.stdout)
        inp = self._fasta()
        proc = self.run_script("blast_sequence_network.py", "-i", inp,
                               "-o", self.tmp / "x", "--identity-threshold", "50")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("fraction", proc.stderr)

    def test_missing_markov_clustering_clean_exit_after_real_blast(self):
        # Everything up to MCL runs for real (blastp 2.16); the missing
        # markov-clustering then aborts cleanly.
        inp = self._fasta()
        prefix = self.tmp / "out"
        proc = self.run_script("blast_sequence_network.py", "-i", inp, "-o", prefix)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("markov-clustering", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        # BLAST output exists; no stray database files remain.
        self.assertTrue(Path(f"{prefix}_blast.txt").exists())
        leftovers = list(self.tmp.glob("*_db.*"))
        self.assertEqual(leftovers, [])

    def test_full_pipeline_with_stubbed_mcl(self):
        # Inject a stub markov_clustering (clustering by connected
        # components) and run main() in-process so the stub is visible;
        # BLAST itself runs for real.
        def fake_run_mcl(matrix, inflation=2.0):
            n = len(matrix)
            parent = list(range(n))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for i in range(n):
                for j in range(n):
                    if matrix[i, j] > 0:
                        parent[find(i)] = find(j)
            groups = {}
            for i in range(n):
                groups.setdefault(find(i), []).append(i)
            return list(groups.values())

        stub = types.ModuleType("markov_clustering")
        stub.run_mcl = fake_run_mcl
        stub.get_clusters = lambda result: result
        had_mc = "markov_clustering" in sys.modules
        saved_mc = sys.modules.get("markov_clustering")
        sys.modules["markov_clustering"] = stub
        inp = self._fasta()
        prefix = self.tmp / "out"
        old_argv = sys.argv
        try:
            sys.argv = [
                "blast_sequence_network.py", "-i", str(inp),
                "-o", str(prefix), "-t", "1",
            ]
            bsn.main()
        finally:
            sys.argv = old_argv
            if had_mc:
                sys.modules["markov_clustering"] = saved_mc
            else:
                sys.modules.pop("markov_clustering", None)

        clusters = Path(f"{prefix}_clusters.txt").read_text()
        self.assertIn("3 members", clusters)   # family A
        self.assertIn("2 members", clusters)   # family B
        self.assertIn("a1, a2, a3", clusters)
        sif = Path(f"{prefix}_network.sif").read_text().splitlines()
        self.assertTrue(sif)
        a_edges = {frozenset(l.split("\t")[::2]) for l in sif}
        self.assertIn(frozenset({"a1", "a2"}), a_edges)
        self.assertIn(frozenset({"b1", "b2"}), a_edges)
        self.assertNotIn(frozenset({"a1", "b1"}), a_edges)
        pdf = Path(f"{prefix}_network.pdf")
        self.assertTrue(pdf.exists() and pdf.stat().st_size > 1000)
        # DB files removed even though MCL was a stub (finally cleanup).
        self.assertEqual(list(self.tmp.glob("*_db.*")), [])

    def test_blast_output_parsing_runs_real_blastp(self):
        inp = self._fasta("tiny.fa")
        prefix = self.tmp / "real"
        # Direct call to run_blast with the module's own resolution.
        out = bsn.run_blast(str(inp), blast_db=str(prefix) + "_db",
                            output_file=str(prefix) + "_blast.txt", evalue=1e-3, threads=1)
        p = Path(out)
        self.assertTrue(p.exists())
        # Self-hits present for all five sequences.
        ids = {line.split("\t")[0] for line in p.read_text().splitlines() if line.strip()}
        self.assertEqual(ids, {"a1", "a2", "a3", "b1", "b2"})


if __name__ == "__main__":
    unittest.main()
