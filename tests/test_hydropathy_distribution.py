"""Functional tests for bin/hydropathy_distribution.py."""

import csv

from common import ScriptTestCase


class HydropathyDistributionTest(ScriptTestCase):
    def read_table(self, out):
        return list(csv.DictReader(out.open()))

    def test_basic_scores(self):
        # "AAAAAA" window 5 -> 2 windows, both mean 1.8
        fa = self.write("p.fa", ">prot\nAAAAAA\n")
        out = self.tmp / "hyd.csv"
        proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out, "-w", "5")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = self.read_table(out)
        self.assertEqual(rows[0]["num_windows"], "2")
        self.assertEqual(rows[0]["mean_score"], "1.8000")
        self.assertEqual(rows[0]["hydropathy_scores"], "1.8000,1.8000")

    def test_known_values_mixed_peptide(self):
        # "AG" window 2: (1.8 + -0.4)/2 = 0.7
        fa = self.write("p.fa", ">prot\nAG\n")
        out = self.tmp / "hyd.csv"
        proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out, "-w", "2")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = self.read_table(out)
        self.assertEqual(rows[0]["mean_score"], "0.7000")

    def test_window_larger_than_sequence_warns_na(self):
        fa = self.write("p.fa", ">short\nAG\n")
        out = self.tmp / "hyd.csv"
        proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out, "-w", "10")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("shorter than window", proc.stderr)
        rows = self.read_table(out)
        self.assertEqual(rows[0]["num_windows"], "0")
        self.assertEqual(rows[0]["mean_score"], "NA")

    def test_invalid_window_clean_error(self):
        fa = self.write("p.fa", ">prot\nAG\n")
        out = self.tmp / "hyd.csv"
        for bad in ("0", "-3"):
            proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out,
                                   "-w", bad)
            self.assertEqual(proc.returncode, 2, f"window={bad}")
            self.assertNotIn("Traceback", proc.stderr, f"window={bad}")

    def test_unknown_residue_warned(self):
        fa = self.write("p.fa", ">prot\nAXAG\n")
        out = self.tmp / "hyd.csv"
        proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out, "-w", "2")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("unknown residue(s) [X]", proc.stderr)

    def test_visualize_pdf_with_slash_in_id(self):
        # regression: '/' in record id crashed savefig
        fa = self.write("p.fa", ">sp/species|prot1\nAGAAGA\n")
        out = self.tmp / "hyd.csv"
        proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out,
                               "-w", "3", "--visualize", "--plot-dir", str(self.tmp))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdf = self.tmp / "hydropathy_sp_species_prot1.pdf"
        self.assertTrue(pdf.exists())
        self.assertTrue(pdf.read_bytes().startswith(b"%PDF"))

    def test_no_plot_for_empty_windows(self):
        fa = self.write("p.fa", ">short\nAG\n")
        out = self.tmp / "hyd.csv"
        proc = self.run_script("hydropathy_distribution.py", "-i", fa, "-o", out,
                               "-w", "10", "--visualize", "--plot-dir", str(self.tmp))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(list(self.tmp.glob("hydropathy_*.pdf")), [])


if __name__ == "__main__":
    unittest.main()
