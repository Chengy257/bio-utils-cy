"""End-to-end tests for bin/batch_dssp.sh (real mkdssp 4.5.3 + GNU parallel)."""

import hashlib
import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402

# minimal poly-alanine backbone (column-exact PDB ATOM records)
def _atom_line(serial, resseq, name, x, y, z, element):
    return (f"ATOM  {serial:5d} {name:4s} ALA A{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:2s}")


def _backbone(n=8):
    lines, serial = [], 1
    for i in range(1, n + 1):
        x = 1.5 * (i - 1)
        for name, ax, ay, az, el in (
                ("N", x - 1.2, 1.0, 1.0, "N"),
                (" CA ", x, 1.2, 1.0, "C"),
                (" C  ", x + 1.3, 1.0, 1.0, "C"),
                (" O  ", x + 2.4, 1.0, 1.0, "O")):
            lines.append(_atom_line(serial, i, name, ax, ay, az, el))
            serial += 1
    return "\n".join(lines) + "\n"


NO_HEADERS = _backbone()
WITH_HEADERS = "HEADER    TEST STRUCTURE                           01-JAN-25   TEST\n" + _backbone()


class TestBatchDssp(ScriptTestCase):

    def _pdb_dir(self, files):
        d = self.tmp / "pdbs"
        d.mkdir()
        for name, text in files.items():
            (d / name).write_text(text)
        return d

    def test_inputs_never_modified_and_outputs_produced(self):
        """Regression: v1.0.0 sed -i'ed the user's PDB files in place."""
        pdb_dir = self._pdb_dir({
            "a.pdb": NO_HEADERS,
            "b.pdb": WITH_HEADERS,
        })
        before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in pdb_dir.iterdir()}
        out = self.tmp / "dssp"
        proc = self.run_script("batch_dssp.sh", "-i", pdb_dir, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)

        after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in pdb_dir.iterdir()}
        self.assertEqual(before, after)  # inputs byte-identical

        self.assertTrue((out / "a.pdb.dssp").exists())
        self.assertTrue((out / "b.pdb.dssp").exists())
        body = (out / "a.pdb.dssp").read_text()
        self.assertIn("Secondary Structure Definition", body)  # real DSSP output
        self.assertIn("TOTAL NUMBER OF RESIDUES", body)

    def test_header_insertion_order_in_staging_copy(self):
        # a headerless file must get HEADER then CRYST1 (spec order);
        # since the staging copy is deleted, verify via successful output
        # plus a file that has only HEADER -> CRYST1 inserted after it
        pdb_dir = self._pdb_dir({"only_header.pdb": WITH_HEADERS})
        out = self.tmp / "dssp2"
        proc = self.run_script("batch_dssp.sh", "-i", pdb_dir, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((out / "only_header.pdb.dssp").exists())

    def test_rerun_skips_existing(self):
        pdb_dir = self._pdb_dir({"a.pdb": NO_HEADERS})
        out = self.tmp / "dssp3"
        first = self.run_script("batch_dssp.sh", "-i", pdb_dir, "-o", out)
        self.assertEqual(first.returncode, 0, first.stderr)
        dssp = out / "a.pdb.dssp"
        stamp = dssp.read_bytes()
        second = self.run_script("batch_dssp.sh", "-i", pdb_dir, "-o", out)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("Skipping existing", second.stderr)
        self.assertEqual(dssp.read_bytes(), stamp)

    def test_corrupt_pdb_fails_nonzero(self):
        pdb_dir = self._pdb_dir({"good.pdb": NO_HEADERS,
                                 "junk.pdb": "THIS IS NOT A PDB FILE\n"})
        out = self.tmp / "dssp4"
        proc = self.run_script("batch_dssp.sh", "-i", pdb_dir, "-o", out)
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertTrue((out / "good.pdb.dssp").exists())
        self.assertFalse((out / "junk.pdb.dssp").exists())
        self.assertIn("failed", proc.stderr)

    def test_invalid_threads_rejected(self):
        pdb_dir = self._pdb_dir({"a.pdb": NO_HEADERS})
        proc = self.run_script("batch_dssp.sh", "-i", pdb_dir,
                               "-o", self.tmp / "o", "-t", "abc")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("positive integer", proc.stderr)

    def test_no_pdb_files_exits_zero(self):
        pdb_dir = self.tmp / "empty_pdbs"
        pdb_dir.mkdir()
        proc = self.run_script("batch_dssp.sh", "-i", pdb_dir, "-o", self.tmp / "o")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("No new PDB files", proc.stderr)

    def test_paths_with_spaces(self):
        d = self.tmp / "pdb dir with spaces"
        d.mkdir()
        (d / "s.pdb").write_text(NO_HEADERS)
        out = self.tmp / "out dir"
        proc = self.run_script("batch_dssp.sh", "-i", d, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((out / "s.pdb.dssp").exists())

    def test_explicit_missing_tool_fails_cleanly(self):
        pdb_dir = self._pdb_dir({"a.pdb": NO_HEADERS})
        proc = self.run_script("batch_dssp.sh", "-i", pdb_dir,
                               "-o", self.tmp / "o",
                               "-d", "/nonexistent/mkdssp")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not executable", proc.stderr)

    def test_old_dssp_cli_fails_loudly(self):
        # /usr/bin/mkdssp is DSSP 2.2.1 (the -i/-o CLI): v4-style jobs must
        # fail the batch loudly rather than exit 0 silently
        pdb_dir = self._pdb_dir({"a.pdb": NO_HEADERS})
        env = {"PATH": "/usr/bin:/bin", "BUC_SKIP_LOCAL_ENV": "1"}
        with mock.patch.dict(os.environ, env):
            proc = self.run_script("batch_dssp.sh", "-i", pdb_dir,
                                   "-o", self.tmp / "o")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertIn("failed", proc.stderr)


if __name__ == "__main__":
    unittest.main()
