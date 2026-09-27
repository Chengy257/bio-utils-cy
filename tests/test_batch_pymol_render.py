"""Tests for bin/batch_pymol_render.py.

The base interpreter lacks the `pymol` module, so the fallback/error paths
run directly; the real render path re-execs under BUC_PYMOL_BIN (pinned in
config/env.local.sh) and is exercised end-to-end.
"""

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402


def _af_style_pdb(n=8):
    """Poly-alanine backbone with varying B-factors (pLDDT)."""
    lines, serial = [], 1
    plddt = [95, 88, 75, 62, 55, 45, 30, 82]
    for i in range(1, n + 1):
        x = 1.5 * (i - 1)
        for name, ax, el in (("N", x - 1.2, "N"), (" CA ", x, "C"),
                             (" C  ", x + 1.3, "C"), (" O  ", x + 2.4, "O")):
            b = plddt[i - 1]
            lines.append(
                f"ATOM  {serial:5d} {name:4s} ALA A{i:4d}    "
                f"{ax:8.3f}   1.200   1.000  1.00 {b:5.2f}           {el:2s}")
            serial += 1
    return "\n".join(lines) + "\nTER\nEND\n"


class TestBatchPymolRender(ScriptTestCase):

    def test_help_without_pymol(self):
        env = dict(os.environ)
        env.pop("BUC_PYMOL_BIN", None)
        with mock.patch.dict(os.environ, env, clear=True):
            proc = self.run_script("batch_pymol_render.py", "--help")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("usage", proc.stdout)

    def test_fallback_error_when_pymol_unavailable(self):
        env = dict(os.environ)
        env.pop("BUC_PYMOL_BIN", None)
        with mock.patch.dict(os.environ, env, clear=True):
            proc = self.run_script("batch_pymol_render.py", "-i", str(self.tmp))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("BUC_PYMOL_BIN", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_no_pdb_files_exits_nonzero(self):
        empty = self.tmp / "no_pdbs"
        empty.mkdir()
        pymol_bin = os.environ.get("BUC_PYMOL_BIN") or str(
            Path.home() / "soft/miniconda3/envs/protein/bin/pymol")
        with mock.patch.dict(os.environ,
                             {"BUC_PYMOL_BIN": pymol_bin,
                              "BUC_SKIP_LOCAL_ENV": "1"}):
            proc = self.run_script("batch_pymol_render.py",
                                   "-i", empty, "-o", self.tmp / "img")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertIn("No PDB files", proc.stderr)

    def test_end_to_end_render_and_exit_code(self):
        pdb_dir = self.tmp / "pdbs"
        pdb_dir.mkdir()
        (pdb_dir / "af1.pdb").write_text(_af_style_pdb())
        out = self.tmp / "img"

        pymol_bin = os.environ.get("BUC_PYMOL_BIN") or str(
            Path.home() / "soft/miniconda3/envs/protein/bin/pymol")
        env = {"BUC_PYMOL_BIN": pymol_bin, "BUC_SKIP_LOCAL_ENV": "1"}
        with mock.patch.dict(os.environ, env):
            proc = self.run_script("batch_pymol_render.py",
                                   "-i", pdb_dir, "-o", out,
                                   "--width", "400", "--height", "300")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        png = out / "af1_plddt.png"
        self.assertTrue(png.exists())
        self.assertGreater(png.stat().st_size, 1000)
        self.assertIn("Rendered 1 / 1", proc.stderr)

    def test_partial_failure_exits_nonzero(self):
        pdb_dir = self.tmp / "pdbs2"
        pdb_dir.mkdir()
        (pdb_dir / "good.pdb").write_text(_af_style_pdb())
        (pdb_dir / "junk.pdb").write_text("NOT A PDB\n")
        out = self.tmp / "img2"

        pymol_bin = os.environ.get("BUC_PYMOL_BIN") or str(
            Path.home() / "soft/miniconda3/envs/protein/bin/pymol")
        with mock.patch.dict(os.environ,
                             {"BUC_PYMOL_BIN": pymol_bin,
                              "BUC_SKIP_LOCAL_ENV": "1"}):
            proc = self.run_script("batch_pymol_render.py",
                                   "-i", pdb_dir, "-o", out,
                                   "--width", "400", "--height", "300")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertTrue((out / "good_plddt.png").exists())
        self.assertFalse((out / "junk_plddt.png").exists())
        self.assertIn("Rendered 1 / 2", proc.stderr)
        self.assertIn("failed", proc.stderr)


if __name__ == "__main__":
    unittest.main()
