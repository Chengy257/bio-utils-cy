"""Offline tests for bin/kegg_enrichment.R.

The CLI surface is testable without network: getopt spec parsing (the
historic -p 0.01 crash), help, required-argument validation, and how
far execution gets before stopping. The terminal point is
environment-dependent — when clusterProfiler is loadable (it needs the
stringi/icu70 LD path) execution proceeds past the package check and
stops at the missing gene-list file; otherwise it stops at the package
check. Both are pinned via a probe.
"""

import os
import subprocess
import unittest

from common import BIN, ScriptTestCase, find_r_pair, find_icu70_lib  # noqa: E402

SCRIPT = "kegg_enrichment.R"

RSCRIPT, R_LIBS = find_r_pair("getopt")
_ICU_LIB = find_icu70_lib()


def _clusterprofiler_loadable():
    if RSCRIPT is None:
        return False
    env = dict(os.environ)
    if R_LIBS:
        env["R_LIBS"] = R_LIBS
    if _ICU_LIB:
        prev = env.get("LD_LIBRARY_PATH")
        env["LD_LIBRARY_PATH"] = _ICU_LIB + (":" + prev if prev else "")
    probe = subprocess.run(
        [RSCRIPT, "-e", "suppressMessages(library(clusterProfiler))"],
        capture_output=True, text=True, timeout=300, env=env,
    )
    return probe.returncode == 0


_CP_LOADABLE = _clusterprofiler_loadable()


@unittest.skipIf(RSCRIPT is None, "R with getopt not available (set BUC_RSCRIPT_BIN)")
class TestKeggEnrichmentCli(ScriptTestCase):

    def _run(self, *args):
        env = dict(os.environ)
        if R_LIBS:
            env["R_LIBS"] = R_LIBS
        if _ICU_LIB:
            prev = env.get("LD_LIBRARY_PATH")
            env["LD_LIBRARY_PATH"] = _ICU_LIB + (":" + prev if prev else "")
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

    def test_value_flags_parse_and_run(self):
        """-p 0.01 used to die inside getopt (flag spec); it must parse."""
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "dosa", "-p", "0.01", "-q", "0.1", "-n", "15", "-a", "BH")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("is not a valid option", proc.stderr)
        if _CP_LOADABLE:
            # Past the package check; stops at the absent gene list.
            self.assertNotIn("Package 'clusterProfiler' is not installed", proc.stderr)
            self.assertIn("Reading gene list", proc.stderr)
        else:
            self.assertIn("Package 'clusterProfiler' is not installed", proc.stderr)

    def test_equals_form_also_parses(self):
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "dosa", "--pvalue=0.01", "--showcat=15")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("is not a valid option", proc.stderr)
        if _CP_LOADABLE:
            self.assertIn("Reading gene list", proc.stderr)
        else:
            self.assertIn("Package 'clusterProfiler' is not installed", proc.stderr)

    def test_libpath_flag_accepted(self):
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "mmu", "-l", "/tmp/does-not-exist-lib")
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("is not a valid option", proc.stderr)
        if _CP_LOADABLE:
            self.assertIn("Reading gene list", proc.stderr)
        else:
            self.assertIn("Package 'clusterProfiler' is not installed", proc.stderr)

    def test_unknown_option_fails_cleanly(self):
        proc = self._run("-i", "genes.txt", "-o", "out", "-t", "map.tsv",
                         "-g", "dosa", "--nope", "1")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
