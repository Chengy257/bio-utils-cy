"""Shared helpers for bin/ script functional tests."""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "bin"


def _rscript_works(rscript):
    """Check that an Rscript binary actually executes code (the wrapper's
    --version can succeed even when the R binary itself is broken)."""
    try:
        proc = subprocess.run([str(rscript), "-e", 'cat("OK")'],
                              capture_output=True, text=True, timeout=120)
        return proc.returncode == 0 and "OK" in proc.stdout
    except Exception:
        return False


def _rscript_has_getopt(rscript):
    try:
        proc = subprocess.run(
            [str(rscript), "-e", 'suppressWarnings(suppressMessages(library(getopt)))'],
            capture_output=True, text=True, timeout=120)
        return proc.returncode == 0
    except Exception:
        return False


def find_rscript():
    """BUC_RSCRIPT_BIN > PATH > conda envs; prefers an R that can load the
    getopt package (project convention for R CLIs), falls back to the first
    runnable R if none has it."""
    candidates = []
    cand = os.environ.get("BUC_RSCRIPT_BIN")
    if cand and Path(cand).exists():
        candidates.append(cand)
    which = shutil.which("Rscript")
    if which:
        candidates.append(which)
    home = Path.home()
    for pattern in ("soft/miniconda3/envs/*/bin/Rscript",
                    "soft/miniconda3/bin/Rscript",
                    "miniconda3/envs/*/bin/Rscript"):
        candidates.extend(sorted(home.glob(pattern)))

    first_working = None
    for c in candidates:
        if not _rscript_works(c):
            continue
        if first_working is None:
            first_working = c
        if _rscript_has_getopt(c):
            return c
    return first_working or which or "Rscript"


RSCRIPT = find_rscript()


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
        if script.endswith(".py"):
            cmd = [sys.executable, str(BIN / script)]
        elif script.endswith(".R"):
            cmd = [RSCRIPT, str(BIN / script)]
        else:
            cmd = ["bash", str(BIN / script)]
        return subprocess.run(cmd + [str(a) for a in args],
                              capture_output=True, text=True, timeout=600)
