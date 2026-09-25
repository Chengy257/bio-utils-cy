"""Functional tests for bin/tissue_specificity_tau.py."""

from common import ScriptTestCase


class TissueSpecificityTauTest(ScriptTestCase):
    def test_headerless_no_pseudo_header(self):
        # regression: headerless input used to emit a bogus "0\ttau" line
        mat = self.write("expr.tsv", "g1\t0\t4\ng2\t2\t2\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = out.read_text().strip().split("\n")
        self.assertEqual(len(lines), 2)  # no header, 2 genes
        self.assertEqual(lines[0].split("\t")[0], "g1")

    def test_tau_values(self):
        # g1 [0,4]: tau=1.0 (single tissue); g2 [2,2]: tau=0.0; g3 [1,2]: tau=0.5
        mat = self.write("expr.tsv", "g1\t0\t4\ng2\t2\t2\ng3\t1\t2\ng4\t0\t0\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = {r.split("\t")[0]: r.split("\t")[1]
                for r in out.read_text().rstrip("\n").split("\n")}
        self.assertEqual(rows["g1"], "1.0")
        self.assertEqual(rows["g2"], "0.0")
        self.assertEqual(rows["g3"], "0.5")
        self.assertEqual(rows["g4"], "")  # all-zero -> NA (empty field)

    def test_with_header_output_carries_header(self):
        mat = self.write("expr.tsv", "gene\tt1\tt2\ng1\t0\t4\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out, "--header")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = out.read_text().strip().split("\n")
        self.assertEqual(lines[0], "gene\ttau")
        self.assertEqual(lines[1].split("\t")[0], "g1")

    def test_negative_values_clean_error(self):
        mat = self.write("expr.tsv", "g1\t-1\t4\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertIn("negative", proc.stderr)

    def test_numeric_gene_ids_verbatim(self):
        mat = self.write("expr.tsv", "10001\t0\t4\n10002\t2\t2\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        first_col = [r.split("\t")[0] for r in out.read_text().strip().split("\n")]
        self.assertEqual(first_col, ["10001", "10002"])

    def test_nan_row_warns(self):
        mat = self.write("expr.tsv", "g1\t0\t4\ng2\t\t2\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("missing values", proc.stderr)

    def test_single_column_error(self):
        mat = self.write("expr.tsv", "g1\ng2\n")
        out = self.tmp / "tau.tsv"
        proc = self.run_script("tissue_specificity_tau.py", "-i", mat, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertNotIn("Traceback", proc.stderr)


if __name__ == "__main__":
    unittest.main()
