"""Tests for bin/plot_structure_3d.py (real mkdssp 4.5.3 available locally)."""

import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import plot_structure_3d as ps3  # noqa: E402
from common import ScriptTestCase  # noqa: E402


def _atom(serial, name, resseq, chain, x, y, z, resn="ALA", element="C", record="ATOM"):
    return (f"{record:6s}{serial:5d} {name:4s} {resn:3s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:2s}")


def _backbone(n=8, chain="A", y_offset=0.0):
    lines, serial = [], 1
    for i in range(1, n + 1):
        x = 1.5 * (i - 1)
        for name, ax, el in (("N", x - 1.2, "N"), (" CA ", x, "C"),
                             (" C  ", x + 1.3, "C"), (" O  ", x + 2.4, "O")):
            lines.append(_atom(serial, name, i, chain, ax, 1.0 + y_offset, 1.0, element=el))
            serial += 1
    return "\n".join(lines) + "\nTER\nEND\n"


class TestParsePdb(ScriptTestCase):

    def test_multichain_same_resnum_kept(self):
        # regression: residues keyed by number only lost chain B's residue 1
        text = (_atom(1, " CA ", 1, "A", 0.0, 0.0, 0.0) + "\n" +
                _atom(2, " CA ", 1, "B", 9.0, 9.0, 9.0) + "\n")
        f = self.write("two.pdb", text + "END\n")
        residues = ps3.parse_pdb(str(f))
        self.assertEqual(len(residues), 2)
        self.assertEqual(sorted(r["chain"] for r in residues), ["A", "B"])

    def test_ca_coordinates_preferred(self):
        text = (_atom(1, " N  ", 1, "A", 0.0, 0.0, 0.0) + "\n" +
                _atom(2, " CA ", 1, "A", 1.5, 1.0, 1.0) + "\n" +
                _atom(3, " C  ", 1, "A", 3.0, 0.0, 0.0) + "\n")
        f = self.write("ca.pdb", text + "END\n")
        r = ps3.parse_pdb(str(f))[0]
        self.assertEqual((r["x"], r["y"], r["z"]), (1.5, 1.0, 1.0))

    def test_hetatm_skipped(self):
        text = (_atom(1, " CA ", 1, "A", 0.0, 0.0, 0.0) + "\n" +
                _atom(2, " O  ", 1, "A", 5.0, 5.0, 5.0, resn="HOH",
                      element="O", record="HETATM") + "\n")
        f = self.write("hoh.pdb", text + "END\n")
        residues = ps3.parse_pdb(str(f))
        self.assertEqual(len(residues), 1)

    def test_sorted_by_chain_then_resnum(self):
        text = (_atom(1, " CA ", 2, "A", 2.0, 0.0, 0.0) + "\n" +
                _atom(2, " CA ", 1, "B", 9.0, 0.0, 0.0) + "\n" +
                _atom(3, " CA ", 1, "A", 1.0, 0.0, 0.0) + "\n")
        f = self.write("srt.pdb", text + "END\n")
        order = [(r["chain"], r["resnum"]) for r in ps3.parse_pdb(str(f))]
        self.assertEqual(order, [("A", 1), ("A", 2), ("B", 1)])


class TestRunDssp(ScriptTestCase):

    def test_ss_map_keyed_by_chain_and_resnum(self):
        f = self.write("bb.pdb", _backbone(8))
        dssp = ps3.detect_dssp()
        if not dssp:
            self.skipTest("mkdssp not available")
        ss_map = ps3.run_dssp(str(f), dssp)
        self.assertIn(("A", 1), ss_map)
        # the poly-ala helix from the backbone fixture is detected on
        # middle residues (regression: chain was not part of the key)
        self.assertEqual(ss_map[("A", 3)], "H")

    def test_multichain_no_key_collision(self):
        # two physically separated chains with continuous numbering; both
        # must appear in the map keyed by (chain, resnum). NOTE: mkdssp
        # 4.5.3's PDB reader mishandles TER records (everything from the
        # second chain is partially dropped), so the fixture omits TER
        # between the chains.
        lines_a = _backbone(8, chain="A").replace("TER\n", "").replace("END\n", "")
        lines_b = []
        serial = 100
        for i in range(101, 109):
            x = 1.5 * (i - 101)
            for name, ax, el in (("N", x - 1.2, "N"), (" CA ", x, "C"),
                                 (" C  ", x + 1.3, "C"), (" O  ", x + 2.4, "O")):
                serial += 1
                lines_b.append(_atom(serial, name, i, "B", ax, 51.0, 1.0, element=el))
        f = self.write("two_chains.pdb", lines_a + "\n".join(lines_b) + "\nEND\n")
        dssp = ps3.detect_dssp()
        if not dssp:
            self.skipTest("mkdssp not available")
        ss_map = ps3.run_dssp(str(f), dssp)
        self.assertIn(("A", 1), ss_map)
        self.assertIn(("B", 101), ss_map)


class TestEndToEnd(ScriptTestCase):

    def test_hydro_mode(self):
        pdb_dir = self.tmp / "pdbs"
        pdb_dir.mkdir()
        (pdb_dir / "p.pdb").write_text(_backbone())
        out = self.tmp / "img"
        proc = self.run_script("plot_structure_3d.py", "-i", pdb_dir, "-o", out,
                               "--color", "hydro", "--width", "400", "--height", "400")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((out / "p.png").exists())
        self.assertGreater((out / "p.png").stat().st_size, 1000)

    def test_ss_mode_with_real_dssp(self):
        pdb_dir = self.tmp / "pdbs2"
        pdb_dir.mkdir()
        (pdb_dir / "p.pdb").write_text(_backbone())
        out = self.tmp / "img2"
        proc = self.run_script("plot_structure_3d.py", "-i", pdb_dir, "-o", out,
                               "--color", "ss", "--width", "400", "--height", "400")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((out / "p.png").exists())

    def test_corrupt_file_fails_batch(self):
        pdb_dir = self.tmp / "pdbs3"
        pdb_dir.mkdir()
        (pdb_dir / "good.pdb").write_text(_backbone())
        (pdb_dir / "junk.pdb").write_text("not a pdb at all\n")
        out = self.tmp / "img3"
        proc = self.run_script("plot_structure_3d.py", "-i", pdb_dir, "-o", out,
                               "--color", "hydro", "--width", "400", "--height", "400")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertTrue((out / "good.png").exists())
        self.assertFalse((out / "junk.png").exists())
        self.assertIn("failed", proc.stderr)

    def test_ss_mode_without_dssp_hard_error(self):
        # in-process: detect_dssp mocked empty; --color ss must hard-error
        pdb_dir = self.tmp / "pdbs4"
        pdb_dir.mkdir()
        (pdb_dir / "p.pdb").write_text(_backbone())
        argv = ["plot_structure_3d.py", "-i", str(pdb_dir),
                "-o", str(self.tmp / "img4"), "--color", "ss"]
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(ps3, "detect_dssp", return_value=""):
            with self.assertLogs(ps3.logger, level="ERROR"):
                with self.assertRaises(SystemExit) as ctx:
                    ps3.main()
        self.assertEqual(ctx.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
