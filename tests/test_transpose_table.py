"""Functional tests for bin/transpose_table.sh."""

from common import ScriptTestCase


class TransposeTableTest(ScriptTestCase):
    def test_basic_transpose(self):
        tsv = self.write("m.tsv", "a\tb\tc\n1\t2\t3\n4\t5\t6\n")
        proc = self.run_script("transpose_table.sh", "-i", tsv)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "a\t1\t4\nb\t2\t5\nc\t3\t6\n")

    def test_double_transpose_roundtrip(self):
        tsv = self.write("m.tsv", "a\tb\tc\n1\t2\t3\n4\t5\t6\n")
        once = self.tmp / "once.tsv"
        twice = self.tmp / "twice.tsv"
        self.assertEqual(self.run_script("transpose_table.sh", "-i", tsv, "-o", once).returncode, 0)
        self.assertEqual(self.run_script("transpose_table.sh", "-i", once, "-o", twice).returncode, 0)
        self.assertEqual(twice.read_text(), tsv.read_text())

    def test_ragged_rows_error(self):
        tsv = self.write("bad.tsv", "a\tb\tc\n1\t2\n")
        proc = self.run_script("transpose_table.sh", "-i", tsv)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("line 2", proc.stderr)
        self.assertIn("expected 3", proc.stderr)

    def test_output_file_atomic(self):
        tsv = self.write("m.tsv", "a\tb\n1\t2\n")
        out = self.tmp / "t.tsv"
        proc = self.run_script("transpose_table.sh", "-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "a\t1\nb\t2\n")
        # no temp leftovers
        leftovers = [p.name for p in self.tmp.iterdir() if p.name.startswith("tmp.")]
        self.assertEqual(leftovers, [])

    def test_failed_run_leaves_no_output(self):
        tsv = self.write("bad.tsv", "a\tb\n1\n")
        out = self.tmp / "t.tsv"
        proc = self.run_script("transpose_table.sh", "-i", tsv, "-o", out)
        self.assertEqual(proc.returncode, 1)
        self.assertFalse(out.exists())

    def test_trailing_empty_field_preserved(self):
        # second row ends with an empty field; streaming must keep it
        tsv = self.write("m.tsv", "a\tb\n1\t\n")
        proc = self.run_script("transpose_table.sh", "-i", tsv)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "a\t1\nb\t\n")

    def test_csv_separator(self):
        csv = self.write("d.csv", "a,b\n1,2\n")
        proc = self.run_script("transpose_table.sh", "-i", csv, "-s", ",")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "a,1\nb,2\n")

    def test_unknown_option_error(self):
        tsv = self.write("m.tsv", "a\n1\n")
        proc = self.run_script("transpose_table.sh", "-i", tsv, "--bogus")
        self.assertEqual(proc.returncode, 1)

    def test_help_and_version(self):
        proc = self.run_script("transpose_table.sh", "-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)
        proc = self.run_script("transpose_table.sh", "-v")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("transpose_table.sh", proc.stdout)


if __name__ == "__main__":
    unittest.main()
