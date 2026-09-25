"""Shared helpers for bin/ script functional tests."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "bin"


class ScriptTestCase(unittest.TestCase):
    """Base class: temp working dir + helpers to run a bin/ script."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, text):
        """Write a (fixture) file into the temp dir and return its path."""
        path = self.tmp / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def run_script(self, script, *args):
        """Run bin/<script> with args; returns CompletedProcess."""
        cmd = [sys.executable if script.endswith(".py") else "bash",
               str(BIN / script)] + [str(a) for a in args]
        if script.endswith(".R"):
            cmd = ["Rscript", str(BIN / script)] + [str(a) for a in args]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=300)
