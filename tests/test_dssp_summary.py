"""Functional tests for bin/dssp_summary.py (synthetic mini .dssp files)."""

import csv

from common import ScriptTestCase

HEADER = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N"


def residue_line(num, aa, ss, acc):
    """Build a fixed-column DSSP residue line: AA@13, SS@16, ACC@34:38."""
    return " " * 13 + aa + "  " + ss + " " * 17 + f"{acc:>4}"


def make_dssp(path):
    lines = [HEADER]
    #          num aa  ss    acc   -> known-AA relative ASA
    lines.append(residue_line(1, "A", "H", 106))    # 106/106 = 1.00
    lines.append(residue_line(2, "R", "H", 62))     # 62/248  = 0.25
    lines.append(residue_line(3, "A", "E", 53))     # 53/106  = 0.50
    lines.append(residue_line(4, "C", "E", 135))    # 135/135 = 1.00
    lines.append(residue_line(5, "X", " ", 50))     # unknown AA: counted, no ASA
    lines.append(residue_line(6, "A", "-", 0))      # 0/106  = 0.00
    path.write_text("\n".join(lines) + "\n")


class DsspSummaryTest(ScriptTestCase):
    def run_summary(self, extra=()):
        out = self.tmp / "summary.csv"
        proc = self.run_script("dssp_summary.py", "-i", self.tmp, "-o", out, *extra)
        return proc, out

    def test_summary_counts_and_ratios(self):
        make_dssp(self.tmp / "prot.dssp")
        proc, out = self.run_summary()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = list(csv.DictReader(out.open()))
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["Filename"], "prot.dssp")
        self.assertEqual(r["TotalResidues"], "6")
        self.assertEqual(r["Helix_Count"], "2")
        self.assertEqual(r["Sheet_Count"], "2")
        self.assertEqual(r["Coil_Count"], "2")     # space + '-' both -> Coil
        self.assertAlmostEqual(float(r["Helix_Ratio"]), 2 / 6, places=3)
        # regression: Raw_space_Count was always 0
        self.assertEqual(r["Raw_space_Count"], "1")
        self.assertEqual(r["Raw_-_Count"], "1")
        self.assertEqual(r["Raw_H_Count"], "2")
        # ASA average over the 5 known-AA residues: (1+.25+.5+1+0)/5
        self.assertAlmostEqual(float(r["AvgRelativeASA"]), 0.55, places=3)

    def test_unknown_aa_counts_towards_ratios(self):
        # regression: unknown-AA residues were dropped from every statistic
        make_dssp(self.tmp / "prot.dssp")
        proc, out = self.run_summary()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        r = list(csv.DictReader(out.open()))[0]
        # Coil includes the X residue (space) -> coil ratio 2/6, sums to 1
        self.assertAlmostEqual(float(r["Coil_Ratio"]), 2 / 6, places=3)
        total = sum(float(r[f"{k}_Ratio"]) for k in ("Helix", "Sheet", "Turn", "Coil"))
        self.assertAlmostEqual(total, 1.0, places=3)

    def test_multiple_files(self):
        make_dssp(self.tmp / "a.dssp")
        make_dssp(self.tmp / "b.dssp")
        proc, out = self.run_summary()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = list(csv.DictReader(out.open()))
        self.assertEqual({r["Filename"] for r in rows}, {"a.dssp", "b.dssp"})

    def test_no_dssp_files_exit_1(self):
        proc, _ = self.run_summary()
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No DSSP files", proc.stderr)

    def test_bad_header_file_skipped(self):
        (self.tmp / "bad.dssp").write_text("not a dssp file\n")
        proc, _ = self.run_summary()
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Could not find data header", proc.stderr)

    def test_chain_breaks_ignored(self):
        make_dssp(self.tmp / "prot.dssp")
        with (self.tmp / "prot.dssp").open("a") as fh:
            fh.write(residue_line(7, "!", "H", 99) + "\n")   # chain break
            fh.write("short\n")                               # too short
        proc, out = self.run_summary()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        r = list(csv.DictReader(out.open()))[0]
        self.assertEqual(r["TotalResidues"], "6")
        self.assertEqual(r["Helix_Count"], "2")


if __name__ == "__main__":
    unittest.main()
