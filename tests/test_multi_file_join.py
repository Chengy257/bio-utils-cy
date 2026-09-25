"""Functional tests for bin/multi_file_join.py."""

from common import ScriptTestCase


class MultiFileJoinTest(ScriptTestCase):
    def _list(self, *paths):
        return self.write("list.txt", "".join(f"{p}\n" for p in paths))

    def test_outer_join_with_headers(self):
        a = self.write("a.tsv", "ID\tvalA\ng1\t0.1\ng2\t0.2\n")
        b = self.write("b.tsv", "ID\tvalB\ng2\t0.3\ng3\t0.4\n")
        lst = self._list(a, b)
        proc, out = self.run_join(lst)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = out.read_text().strip().split("\n")
        self.assertEqual(rows[0], "ID\ta.tsv\tb.tsv")
        self.assertEqual(len(rows), 4)  # header + g1, g2, g3 (outer join)
        by_id = {r.split("\t")[0]: r.split("\t")[1:] for r in rows[1:]}
        self.assertEqual(by_id["g1"], ["0.1", "NA"])
        self.assertEqual(by_id["g2"], ["0.2", "0.3"])
        self.assertEqual(by_id["g3"], ["NA", "0.4"])

    def run_join(self, file_list, extra=(), columns="1,2"):
        out = self.tmp / "merged.tsv"
        proc = self.run_script("multi_file_join.py", "-i", file_list,
                               "-c", columns, "-o", out, *extra)
        return proc, out

    def test_no_header_keeps_first_data_row(self):
        a = self.write("a.tsv", "g1\t10\ng2\t20\n")
        b = self.write("b.tsv", "g1\tx\ng2\ty\n")
        lst = self._list(a, b)
        proc, out = self.run_join(lst, extra=["--no-header"])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = out.read_text().strip().split("\n")
        self.assertEqual(len(rows), 3)  # header + g1 + g2: first data row kept
        by_id = {r.split("\t")[0]: r.split("\t")[1:] for r in rows[1:]}
        self.assertEqual(by_id["g1"], ["10", "x"])

    def test_numeric_ids_join_across_type_inference(self):
        # a has pure-numeric IDs (pandas would infer int64), b has mixed
        # (object). Values must still join by string identity.
        a = self.write("a.tsv", "1\t10\n2\t20\n")
        b = self.write("b.tsv", "1\tx\n2\ty\nQ\tz\n")
        lst = self._list(a, b)
        proc, out = self.run_join(lst, extra=["--no-header"])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        ids = [r.split("\t")[0] for r in out.read_text().strip().split("\n")[1:]]
        self.assertEqual(sorted(ids), ["1", "2", "Q"])
        rows = out.read_text().strip().split("\n")[1:]
        by_id = {r.split("\t")[0]: r.split("\t")[1:] for r in rows}
        self.assertEqual(by_id["1"], ["10", "x"])

    def test_missing_file_writes_output_but_exits_1(self):
        a = self.write("a.tsv", "ID\tvalA\ng1\t0.1\n")
        ghost = self.tmp / "ghost.tsv"
        lst = self._list(a, ghost)
        proc, out = self.run_join(lst)
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertIn("failed to load", proc.stderr)
        self.assertTrue(out.exists())
        self.assertIn("g1", out.read_text())

    def test_bad_columns_clean_error(self):
        a = self.write("a.tsv", "ID\tv\ng1\t1\n")
        lst = self._list(a)
        for bad in ("1,a", "1", "0,2", "1,2,3"):
            proc, _ = self.run_join(lst, columns=bad)
            self.assertNotEqual(proc.returncode, 0, f"columns={bad!r}")
            self.assertNotIn("Traceback", proc.stderr, f"columns={bad!r}")

    def test_duplicate_basename_gets_unique_column(self):
        a = self.write("x/f.tsv", "ID\tv1\ng1\t1\n")
        b = self.write("y/f.tsv", "ID\tv2\ng1\t2\n")
        lst = self._list(a, b)
        proc, out = self.run_join(lst)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        header = out.read_text().strip().split("\n")[0]
        self.assertIn("f.tsv", header)
        self.assertIn("f_2.tsv", header)

    def test_duplicate_id_warns(self):
        a = self.write("a.tsv", "ID\tv\ng1\t1\ng1\t2\n")
        b = self.write("b.tsv", "ID\tw\ng1\t3\n")
        lst = self._list(a, b)
        proc, out = self.run_join(lst)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("duplicate ID", proc.stderr)
        # cartesian product: 2 g1 rows from a x 1 g1 row from b = 2 data rows
        self.assertEqual(len(out.read_text().strip().split("\n")), 3)

    def test_sep_option(self):
        a = self.write("a.csv", "ID,valA\ng1,0.1\n")
        lst = self._list(a)
        proc, out = self.run_join(lst, extra=["--sep", ","])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text().strip().split("\n")[0], "ID,a.csv")


if __name__ == "__main__":
    unittest.main()
