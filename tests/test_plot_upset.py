"""Functional tests for bin/plot_upset.R (real UpSetR rendering)."""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "UpSetR")


@unittest.skipUnless(RS, "R with getopt+UpSetR not available")
class TestPlotUpset(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # relative output paths land in the temp dir
        try:
            return subprocess.run([RS, str(BIN / "plot_upset.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def _listfile(self, *paths):
        return self.write("lists.txt", "\n".join(str(p) for p in paths) + "\n")

    def test_basic_render(self):
        a = self.write("A.txt", "G1\nG2\nG3\nG4\n")
        b = self.write("B.txt", "G3\nG4\nG5\n")
        lst = self._listfile(a, b)
        out = self.tmp / "res"
        proc = self._run("-i", lst, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdf = self.tmp / "res_upset_plot.pdf"
        self.assertTrue(pdf.exists() and pdf.stat().st_size > 2000, proc.stderr)
        self.assertIn("Number of sets:        2", proc.stderr)
        self.assertIn("Total unique genes:    5", proc.stderr)

    def test_crlf_lists_still_intersect(self):
        # regression: CRLF IDs kept a trailing \r and the intersection
        # against an LF list was silently empty
        a = self.write("A.txt", "G1\r\nG2\r\nG3\r\n")
        b = self.write("B.txt", "G1\nG2\nG4\n")
        lst = self._listfile(a, b)
        proc = self._run("-i", lst, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Total unique genes:    4", proc.stderr)  # G1,G2 common

    def test_cr_only_line_not_a_gene(self):
        a = self.write("A.txt", "G1\r\n\r\nG2\r\n")
        b = self.write("B.txt", "G1\nG2\n")
        lst = self._listfile(a, b)
        proc = self._run("-i", lst, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Read 2 unique genes from", proc.stderr)

    def test_duplicate_basenames_disambiguated(self):
        # regression: dirA/genes.txt and dirB/genes.txt silently merged
        d1 = self.tmp / "dirA"
        d2 = self.tmp / "dirB"
        d1.mkdir()
        d2.mkdir()
        (d1 / "genes.txt").write_text("G1\nG2\n")
        (d2 / "genes.txt").write_text("G3\nG4\n")
        lst = self._listfile(d1 / "genes.txt", d2 / "genes.txt")
        proc = self._run("-i", lst, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("renamed to 'genes_2'", proc.stderr)
        self.assertIn("Number of sets:        2", proc.stderr)
        self.assertIn("genes_2                        2 genes", proc.stderr)

    def test_hidden_intersections_reported(self):
        # 3 sets over 4 genes: patterns {A},{B},{C},{ABC} -> 4 intersections;
        # display only 2
        lists = []
        for i, g in enumerate((["G1", "G4"], ["G2", "G4"], ["G3", "G4"])):
            p = self.tmp / f"S{i}.txt"
            p.write_text("\n".join(g) + "\n")
            lists.append(p)
        lst = self._listfile(*lists)
        proc = self._run("-i", lst, "-o", "res", "-n", "2")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("showing the 2 largest", proc.stderr)

    def test_empty_list_skipped(self):
        a = self.write("A.txt", "G1\nG2\n")
        e = self.write("E.txt", "\n")
        b = self.write("B.txt", "G2\nG3\n")
        lst = self._listfile(a, e, b)
        proc = self._run("-i", lst, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("is empty. Skipping", proc.stderr)

    def test_single_nonempty_list_fails(self):
        a = self.write("A.txt", "G1\nG2\n")
        e = self.write("E.txt", "")
        lst = self._listfile(a, e)
        proc = self._run("-i", lst, "-o", "res")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("at least 2 non-empty gene lists", proc.stderr)

    def test_missing_list_file_fails(self):
        a = self.write("A.txt", "G1\n")
        b = self.write("B.txt", "G1\n")
        lst = self._listfile(a, self.tmp / "nope.txt")
        proc = self._run("-i", lst, "-o", "res")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("were not found", proc.stderr)

    def test_bad_order_by_rejected(self):
        a = self.write("A.txt", "G1\n")
        b = self.write("B.txt", "G1\n")
        lst = self._listfile(a, b)
        proc = self._run("-i", lst, "-o", "res", "--order-by", "size")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("'degree' or 'freq'", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
