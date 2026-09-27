"""Functional tests for bin/plot_go_heatmap_batch.sh.

The bash layer is exercised through a stub R script that honours the
-f/-g/-o contract of plot_heatmap_multi.R and dumps everything it
received into STUB_DUMP/call.txt. One integration test runs the real
plot_heatmap_multi.R end to end (skipped without ComplexHeatmap).
"""

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

STUB_R = r"""
args <- commandArgs(trailingOnly = TRUE)
get1 <- function(flag) {
    i <- which(args == flag)
    if (length(i) == 0) NA_character_ else args[i + 1]
}
fl <- get1("-f"); gr <- get1("-g"); op <- get1("-o")
if (is.na(fl) || is.na(gr) || is.na(op)) {
    cat("STUB: expected -f/-g/-o, got:", paste(args, collapse = " "), "\n")
    quit(status = 2)
}
dump <- Sys.getenv("STUB_DUMP")
lines <- c(paste0("n_files=", length(readLines(fl))))
for (p in readLines(fl)) {
    lines <- c(lines, paste0(basename(p), "=", length(readLines(p))))
}
lines <- c(lines, paste0("groups=", paste(readLines(gr), collapse = ",")))
lines <- c(lines, paste0("prefix=", op))
writeLines(lines, file.path(dump, "call.txt"))
"""


class TestPlotGoHeatmapBatch(ScriptTestCase):
    def setUp(self):
        super().setUp()
        self.dump = self.tmp / "dump"
        self.dump.mkdir()
        self.r_stub = self.write("stub.R", STUB_R)

    def _batch(self, *args):
        env = dict(os.environ, STUB_DUMP=str(self.dump))
        with mock.patch.dict(os.environ, env):
            return self.run_script("plot_go_heatmap_batch.sh", "-r", self.r_stub, *args)

    def _call(self):
        return (self.dump / "call.txt").read_text().splitlines()

    GMT = ("GO:1\tmRNA.splicing\tG1\tG2\tG3\n"
           "GO:2\tregulation\tG4\tG5\t\n")  # trailing tab -> empty field

    def _fixture(self, gmt_text=None, columns=("S3", "S1", "S4", "S2")):
        """8-gene matrix with columns deliberately ordered differently
        from the sample-info rows (S1=Control, S2=Treat)."""
        gmt = self.write("go.gmt", self.GMT if gmt_text is None else gmt_text)
        header = "gene\t" + "\t".join(columns)
        rows = "\n".join(f"G{i}\t{i}1\t{i}2\t{i}3\t{i}4" for i in range(1, 9))
        expr = self.write("expr.tsv", header + "\n" + rows + "\n")
        info = self.write("info.csv",
                          "SampleID,Group\nS1,Control\nS2,Treat\n"
                          "S3,Control\nS4,Treat\n")
        ids = self.write("ids.txt", "GO:1\nGO:2\n")
        return gmt, expr, info, ids

    def test_groups_follow_matrix_column_order(self):
        gmt, expr, info, ids = self._fixture()
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        groups = [ln for ln in self._call() if ln.startswith("groups=")][0]
        # matrix columns S3,S1,S4,S2 -> Control,Control,Treat,Treat
        self.assertEqual(groups, "groups=Control,Control,Treat,Treat")

    def test_subset_rows_exact_and_trailing_tab_ignored(self):
        gmt, expr, info, ids = self._fixture()
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        counts = {ln.split("=")[0]: ln.split("=")[1] for ln in self._call() if "=" in ln}
        self.assertEqual(counts.get("n_files"), "2")
        # header + exactly the GO genes; the trailing-tab empty field of
        # GO:2 must NOT widen the match to the whole matrix
        self.assertEqual(counts.get("tempExp.GO_1_mRNA_splicing"), "4")
        self.assertEqual(counts.get("tempExp.GO_2_regulation"), "3")

    def test_dot_in_description_sanitized_and_unique(self):
        gmt, expr, info, ids = self._fixture(
            "GO:1\tmRNA.splicing\tG1\n"
            "GO:2\tmRNA.processing\tG2\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        names = [ln.split("=")[0] for ln in self._call() if ln.startswith("tempExp.")]
        self.assertEqual(sorted(names),
                         ["tempExp.GO_1_mRNA_splicing", "tempExp.GO_2_mRNA_processing"])

    def test_regex_dot_does_not_overmatch(self):
        gmt, expr, info, ids = self._fixture("GO:1\ttest\tAT1G01010.1\n")
        expr = self.write("expr.tsv",
                          "gene\tS1\tS2\n"
                          "AT1G01010.1\t5\t6\n"
                          "AT1G01010X1\t7\t8\n"
                          "G3\t1\t2\n")
        info = self.write("info.csv", "SampleID,Group\nS1,A\nS2,B\n")
        ids = self.write("ids.txt", "GO:1\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("1/1 genes matched", proc.stderr)
        self.assertIn("tempExp.GO_1_test=2", self._call())

    def test_word_boundary_neighbour_not_matched(self):
        gmt, expr, info, ids = self._fixture("GO:1\tkinase\tbHLH38\n")
        expr = self.write("expr.tsv",
                          "gene\tS1\tS2\n"
                          "bHLH38\t5\t6\n"
                          "bHLH38-like\t7\t8\n")
        info = self.write("info.csv", "SampleID,Group\nS1,A\nS2,B\n")
        ids = self.write("ids.txt", "GO:1\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("tempExp.GO_1_kinase=2", self._call())

    def test_header_line_never_rematched(self):
        # a gene literally named like a sample must not duplicate the
        # header row into the subset (old code matched the header too)
        gmt, expr, info, ids = self._fixture("GO:1\tsample\tS1\n")
        expr = self.write("expr.tsv", "gene\tS1\tS2\nS1\t5\t6\nG2\t7\t8\n")
        info = self.write("info.csv", "SampleID,Group\nS1,A\nS2,B\n")
        ids = self.write("ids.txt", "GO:1\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("tempExp.GO_1_sample=2", self._call())

    def test_crlf_inputs_processed(self):
        gmt, expr, info, ids = self._fixture("GO:1\tcrlf.term\tG1\tG2\r\n",
                                             columns=("S1", "S2"))
        info = self.write("info_crlf.csv", "SampleID,Group\r\nS1,A\r\nS2,B\r\n")
        ids = self.write("ids.txt", "GO:1\r\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("2/2 genes matched", proc.stderr)

    def test_duplicate_go_ids_processed_once(self):
        gmt, expr, info, ids = self._fixture("GO:1\tonce\tG1\n", columns=("S1", "S2"))
        ids = self.write("ids.txt", "GO:1\nGO:1\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("n_files=1", self._call())

    def test_sample_without_group_fails(self):
        gmt, expr, info, ids = self._fixture()
        info = self.write("bad.csv",
                          "SampleID,Group\nS1,Control\nS2,Treat\nS3,Control\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("S4 has no group", proc.stderr)

    def test_go_missing_from_gmt_skips(self):
        gmt, expr, info, ids = self._fixture()
        ids = self.write("ids.txt", "GO:1\nGO:9999\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("GO:9999' not found in GMT", proc.stderr)
        self.assertIn("Processed 1 GO terms (1 skipped)", proc.stderr)

    def test_all_go_skipped_exits_1(self):
        gmt, expr, info, ids = self._fixture()
        ids = self.write("ids.txt", "GO:9998\nGO:9999\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No GO terms produced valid gene subsets", proc.stderr)

    def test_no_gene_matches_matrix_skips_term(self):
        gmt, expr, info, ids = self._fixture("GO:1\tnope\tG99\tG98\n")
        ids = self.write("ids.txt", "GO:1\n")
        proc = self._batch("-g", gmt, "-o", "out/hm", ids, expr, info)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No genes from GO:1", proc.stderr)

    def test_cli_errors(self):
        gmt, expr, info, ids = self._fixture()
        proc = self._batch("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)

        proc = self._batch(ids, expr, info)  # missing -g
        self.assertEqual(proc.returncode, 1)
        self.assertIn("-g GMT_FILE is required", proc.stderr)

        proc = self._batch("-g", gmt, ids, expr)  # 2 positionals
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Expected exactly 3 positional", proc.stderr)

        proc = self._batch("-g", self.tmp / "nope.gmt", "-o", "out/hm",
                           ids, expr, info)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Required file not found", proc.stderr)

    def test_full_chain_real_r(self):
        rs, rlibs = find_r_pair("getopt", "ComplexHeatmap", "circlize")
        if not rs:
            self.skipTest("R with getopt+ComplexHeatmap not available")
        gmt, expr, info, ids = self._fixture()
        outdir = self.tmp / "out"
        outdir.mkdir()
        env = dict(os.environ)
        env["BUC_RSCRIPT_BIN"] = rs
        if rlibs:
            env["R_LIBS"] = rlibs
        with mock.patch.dict(os.environ, env):
            proc = self.run_script("plot_go_heatmap_batch.sh",
                                   "-g", gmt, "-o", str(outdir / "hm"),
                                   ids, expr, info)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        pdfs = sorted(outdir.glob("*.pdf"))
        self.assertEqual(len(pdfs), 3)  # 2 single + 1 multi
        for p in pdfs:
            self.assertGreater(p.stat().st_size, 3000, p)


if __name__ == "__main__":
    unittest.main()
