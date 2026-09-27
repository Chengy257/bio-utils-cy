"""Offline tests for bin/kegg_enrichment.R.

clusterProfiler is not installed on this machine, so the enrichment run
itself cannot execute — but everything up to the package check is real
and testable: getopt spec parsing (the -p 0.01 crash), help, required-
argument validation, and the package check message.
"""

import os
import subprocess
import unittest

from common import BIN, ScriptTestCase, find_r_pair  # noqa: E402

SCRIPT = "kegg_enrichment.R"

RSCRIPT, R_LIBS = find_r_pair("getopt")


@unittest.skipIf(RSCRIPT is None, "R with getopt not available (set BUC_RSCRIPT_BIN)")
class TestKeggEnrichmentCli(ScriptTestCase):

    def _run(self, *args):
        env = dict(os.environ)
        if R_LIBS:
            env["R_LIBS"] = R_LIBS
        return subprocess.run([RSCRIPT, str(BIN / SCRIPT)] + [str(a) for a in args],
                              capture_output=True, text=True, timeout=300, env=env)

    def test_help_on_stdout(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Usage:", proc.stdout)
        self.assertIn("KEGG organism code", proc.stdout)
        self.assertIn("-p 0.01", proc.stdout)  # usage example
        self.assertEqual(proc.stderr, "")

    def test_missing_required_args_fail(self):
        proc = self._run()
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Missing required arguments", proc.stderr)

    def test_value_flags_parse_and_reach_package_check(self):
        """-p 0.01 used to die inside getopt (flag spec); it must parse."""
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "dosa", "-p", "0.01", "-q", "0.1", "-n", "15", "-a", "BH")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("is not a valid option", proc.stderr)
        self.assertIn("Package 'clusterProfiler' is not installed", proc.stderr)

    def test_equals_form_also_parses(self):
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "dosa", "--pvalue=0.01", "--showcat=15")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("is not a valid option", proc.stderr)
        self.assertIn("Package 'clusterProfiler' is not installed", proc.stderr)

    def test_libpath_flag_accepted(self):
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "mmu", "-l", "/tmp/does-not-exist-lib")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("is not a valid option", proc.stderr)
        self.assertIn("Package 'clusterProfiler' is not installed", proc.stderr)

    def test_unknown_option_fails_cleanly(self):
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "dosa", "--nope", "1")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
