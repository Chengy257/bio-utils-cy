"""End-to-end tests for bin/deg_wilcox_test.R (edgeR TMM + Wilcoxon).

Runs the real script with the R found by find_r_pair("edgeR", "getopt")
and synthetic 10-gene matrices. Fully offline.
"""

import os
import subprocess
import unittest
from pathlib import Path

from common import BIN, ScriptTestCase, find_r_pair  # noqa: E402

SCRIPT = "deg_wilcox_test.R"
NEEDED_PKGS = ("edgeR", "getopt")

RSCRIPT, R_LIBS = find_r_pair(*NEEDED_PKGS)


def make_counts(n_tied=0):
    """6 samples (a1..a3, b1..b3). g_up is higher in b, g_down lower;
    background genes vary; optional fully-tied genes (NaN p-values)."""
    rows = ["gene\ta1\ta2\ta3\tb1\tb2\tb3"]
    rows.append("g_up\t100\t120\t110\t1000\t1200\t1100")
    rows.append("g_down\t900\t880\t910\t100\t90\t110")
    rows.append("g_bg1\t50\t60\t55\t52\t58\t54")
    rows.append("g_bg2\t30\t28\t32\t31\t29\t33")
    rows.append("g_low\t2\t1\t2\t1\t2\t1")  # rowSums=9 -> filtered by default
    for i in range(n_tied):
        rows.append(f"g_tied{i}\t7\t7\t7\t7\t7\t7")
    return "\n".join(rows) + "\n"


CONDITIONS = "a\ta\ta\tb\tb\tb\n"


@unittest.skipIf(RSCRIPT is None,
                 f"R with {', '.join(NEEDED_PKGS)} not available (set BUC_RSCRIPT_BIN)")
class TestDegWilcoxTest(ScriptTestCase):

    def _env(self):
        env = dict(os.environ)
        if R_LIBS:
            env["R_LIBS"] = R_LIBS
        return env

    def _run(self, *args):
        return subprocess.run([RSCRIPT, str(BIN / SCRIPT)] + [str(a) for a in args],
                              capture_output=True, text=True, timeout=600,
                              env=self._env())

    def _fixtures(self, counts=None, conditions=CONDITIONS):
        cpath = self.write("counts.tsv", counts or make_counts())
        tpath = self.write("conditions.tsv", conditions)
        return cpath, tpath

    def _read_tsv(self, path):
        lines = Path(path).read_text().splitlines()
        header = lines[0].split("\t")
        return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

    def test_end_to_end_direction_and_fdr(self):
        cpath, tpath = self._fixtures()
        prefix = self.tmp / "out"
        proc = self._run("-c", cpath, "-t", tpath, "-o", prefix)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "")  # logs on stderr, help only on stdout

        header, rows = self._read_tsv(str(prefix) + "_all.tsv")
        self.assertEqual(header, ["gene_id", "log2FoldChange", "pvalue", "padj"])
        by_gene = {r["gene_id"]: r for r in rows}
        self.assertEqual(len(rows), 4)  # g_low filtered (rowSums=9 <= 10)
        self.assertNotIn("g_low", by_gene)

        # default reference = first condition in file order (a): up in b -> FC > 0
        self.assertGreater(float(by_gene["g_up"]["log2FoldChange"]), 1.0)
        self.assertLess(float(by_gene["g_down"]["log2FoldChange"]), -1.0)
        # 3 vs 3 Wilcoxon cannot get below p~0.08; p-values must be finite
        self.assertLess(float(by_gene["g_up"]["pvalue"]), 0.1)
        self.assertTrue(all(r["padj"] not in ("NA", "NaN", "") for r in rows))

        # padj equals BH over the reported p-values
        ps = [float(r["pvalue"]) for r in rows]
        ranks = sorted(range(len(ps)), key=lambda i: ps[i])
        bh = [None] * len(ps)
        prev = 1.0
        for rank_pos, i in enumerate(reversed(ranks), 1):
            k = len(ps) - rank_pos + 1
            val = min(prev, ps[i] * len(ps) / k)
            bh[i] = val
            prev = val
        for r, expected in zip(rows, bh):
            self.assertAlmostEqual(float(r["padj"]), round(expected, 6), places=5)

    def test_degs_file_with_loose_fdr(self):
        """With 3v3 samples no gene reaches padj<0.05; -p 0.5 exercises the
        _DEGs.tsv mechanics (g_up/g_down in, background out)."""
        cpath, tpath = self._fixtures()
        prefix = self.tmp / "out"
        proc = self._run("-c", cpath, "-t", tpath, "-o", prefix, "-p", "0.5")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, sig = self._read_tsv(str(prefix) + "_DEGs.tsv")
        sig_ids = {r["gene_id"] for r in sig}
        self.assertIn("g_up", sig_ids)
        self.assertIn("g_down", sig_ids)
        self.assertNotIn("g_bg1", sig_ids)
        self.assertIn("Significant (FDR <0.5): 2", proc.stderr)

    def test_control_flag_flips_direction(self):
        cpath, tpath = self._fixtures()
        prefix = self.tmp / "out"
        proc = self._run("-c", cpath, "-t", tpath, "-o", prefix, "-r", "b")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        _, rows = self._read_tsv(str(prefix) + "_all.tsv")
        by_gene = {r["gene_id"]: r for r in rows}
        # reference is now b: g_up is DOWN relative to b
        self.assertLess(float(by_gene["g_up"]["log2FoldChange"]), -1.0)
        self.assertGreater(float(by_gene["g_down"]["log2FoldChange"]), 1.0)
        self.assertIn("vs b", proc.stderr.replace("(n=3)", "").replace(" ", " "))
        self.assertIn("[INFO] Comparing: a", proc.stderr)

    def test_tied_genes_reported_not_silent(self):
        """Fully tied CPM values give p=NaN (verified: raw counts only tie
        BEFORE normalization — TMM breaks raw ties unless every library has
        the same composition). All-equal columns -> TMM factors exactly 1.
        """
        tied_rows = ["gene\ta1\ta2\ta3\tb1\tb2\tb3",
                     "g_tied0\t100\t100\t100\t100\t100\t100",
                     "g_tied1\t50\t50\t50\t50\t50\t50",
                     "g_bg1\t425\t425\t425\t425\t425\t425",
                     "g_bg2\t425\t425\t425\t425\t425\t425"]
        cpath = self.write("counts.tsv", "\n".join(tied_rows) + "\n")
        tpath = self.write("conditions.tsv", CONDITIONS)
        prefix = self.tmp / "out"
        proc = self._run("-c", cpath, "-t", tpath, "-o", prefix)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("4 gene(s) have undefined p-values", proc.stderr)

        _, rows = self._read_tsv(str(prefix) + "_all.tsv")
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(r["padj"] in ("NA", "NaN", "") for r in rows))

        _, sig = self._read_tsv(str(prefix) + "_DEGs.tsv")
        self.assertEqual(len(sig), 0)  # NA padj never passes the FDR filter

    def test_min_count_parameter(self):
        cpath, tpath = self._fixtures()
        prefix = self.tmp / "o1"
        self._run("-c", cpath, "-t", tpath, "-o", prefix)
        _, rows = self._read_tsv(str(prefix) + "_all.tsv")
        self.assertNotIn("g_low", {r["gene_id"] for r in rows})  # rowSums=9

        prefix = self.tmp / "o2"
        proc = self._run("-c", cpath, "-t", tpath, "-o", prefix, "-m", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("rowSums > 1", proc.stderr)
        _, rows = self._read_tsv(str(prefix) + "_all.tsv")
        self.assertIn("g_low", {r["gene_id"] for r in rows})

    def test_numeric_condition_labels_stay_character(self):
        cpath = self.write("counts.tsv", make_counts())
        tpath = self.write("conditions.tsv", "1\t1\t1\t2\t2\t2\n")
        prefix = self.tmp / "out"
        proc = self._run("-c", cpath, "-t", tpath, "-o", prefix)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("[INFO] Comparing: 2", proc.stderr)
        _, rows = self._read_tsv(str(prefix) + "_all.tsv")
        by_gene = {r["gene_id"]: r for r in rows}
        self.assertGreater(float(by_gene["g_up"]["log2FoldChange"]), 1.0)

    def test_condition_count_mismatch_fails(self):
        cpath = self.write("counts.tsv", make_counts())
        tpath = self.write("conditions.tsv", "a\ta\tb\tb\tb\tb\tb\n")
        proc = self._run("-c", cpath, "-t", tpath, "-o", self.tmp / "out")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("condition count", proc.stderr)

    def test_unknown_control_fails(self):
        cpath, tpath = self._fixtures()
        proc = self._run("-c", cpath, "-t", tpath, "-o", self.tmp / "out", "-r", "nosuch")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not found", proc.stderr)
        self.assertIn("a, b", proc.stderr)

    def test_missing_files_fail(self):
        cpath = self.write("counts.tsv", make_counts())
        proc = self._run("-c", cpath, "-t", self.tmp / "nope.tsv", "-o", self.tmp / "out")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Condition file not found", proc.stderr)

    def test_help_on_stdout(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Usage: Rscript deg_wilcox_test.R", proc.stdout)
        self.assertIn("--control", proc.stdout)
        self.assertIn("--min_count", proc.stdout)


if __name__ == "__main__":
    unittest.main()
