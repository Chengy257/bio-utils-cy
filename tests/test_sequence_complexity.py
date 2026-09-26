"""Functional tests for bin/sequence_complexity.py.

The seq_entropies package is not installed in this environment, so
coverage covers validation paths and the clean dependency abort; the
entropy computation itself stays untested here (same as pre-existing
smoke coverage).
"""

import unittest

from common import ScriptTestCase


class SequenceComplexityTest(ScriptTestCase):

    def _fasta(self, text, name="in.fa"):
        return self.write(name, text)

    def test_missing_dependency_clean_exit(self):
        inp = self._fasta(">s1\nACGTACGTAC\n")
        proc = self.run_script("sequence_complexity.py", "-i", inp, "-o", self.tmp / "o.tsv")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("seq_entropies", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_help_and_version_without_dependency(self):
        proc = self.run_script("sequence_complexity.py", "--help")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("nats", proc.stdout)
        proc = self.run_script("sequence_complexity.py", "--version")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("1.1.0", proc.stdout)

    def test_window_validation(self):
        inp = self._fasta(">s1\nACGTACGTAC\n")
        proc = self.run_script("sequence_complexity.py", "-i", inp,
                               "-o", self.tmp / "o.tsv", "-w", "0")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Window size", proc.stderr)

    def test_base_validation(self):
        inp = self._fasta(">s1\nACGTACGTAC\n")
        for bad in ("0", "-2", "1"):
            proc = self.run_script("sequence_complexity.py", "-i", inp,
                                   "-o", self.tmp / "o.tsv", "-b", bad)
            self.assertEqual(proc.returncode, 1, bad)
            self.assertIn("Base", proc.stderr)

    def test_base_accepts_float(self):
        # e as a float literal must pass argparse (previously type=int).
        inp = self._fasta(">s1\nACGTACGTAC\n")
        proc = self.run_script("sequence_complexity.py", "-i", inp,
                               "-o", self.tmp / "o.tsv", "-b", "2.718281828")
        # Fails later on the missing dependency, not at argparse.
        self.assertEqual(proc.returncode, 1)
        self.assertIn("seq_entropies", proc.stderr)

    def test_missing_input(self):
        proc = self.run_script("sequence_complexity.py", "-i", self.tmp / "no.fa",
                               "-o", self.tmp / "o.tsv")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not found", proc.stderr)


if __name__ == "__main__":
    unittest.main()
