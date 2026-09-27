"""End-to-end tests for bin/merge_featurecounts.py.

All-Python, fully offline: hand-written .count fixtures exercising the
merge contract, the outer-join NaN policy (counts -> 0, FPKM/TPM -> NaN +
warning), the effLength audit, and the *.fc.tsv rejection.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402


def count_tsv(genes, counts=None, fpkm=None, tpm=None, eff=None):
    """Build a .count table: genes = [id, ...]; per-gene overrides by dict."""
    rows = ["id\teffLength\tcounts\tfpkm\ttpm"]
    for i, gid in enumerate(genes):
        e = (eff or {}).get(gid, 100 + 10 * i)
        c = (counts or {}).get(gid, 10 + i)
        f = (fpkm or {}).get(gid, round(1.5 + 0.1 * i, 3))
        t = (tpm or {}).get(gid, round(2.5 + 0.1 * i, 3))
        rows.append(f"{gid}\t{e}\t{c}\t{f}\t{t}")
    return "\n".join(rows) + "\n"


class TestMergeFeaturecounts(ScriptTestCase):

    def _run(self, *args):
        return self.run_script("merge_featurecounts.py", *args)

    def test_basic_merge_three_matrices_and_logs(self):
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.count").write_text(count_tsv(["g1", "g2"]))
        (indir / "S2.count").write_text(count_tsv(["g1", "g2"]))
        (indir / "S1.log").write_text("Assigned\t100\nUnassigned_Ambiguity\t2\n")

        proc = self._run("-i", str(indir), "-o", str(self.tmp / "out"))
        self.assertEqual(proc.returncode, 0, proc.stderr)

        out = self.tmp / "out"
        counts = (out / "count.matrix.tsv").read_text().splitlines()
        self.assertEqual(counts[0], "id\tS1\tS2")
        self.assertEqual(counts[1], "g1\t10\t10")
        self.assertEqual(counts[2], "g2\t11\t11")

        fpkm = (out / "GeneExpression_FPKM.xls").read_text().splitlines()
        self.assertEqual(fpkm[1].split("\t")[1:], ["1.5", "1.5"])

        efflen = (out / "effLength.txt").read_text().splitlines()
        self.assertEqual(efflen, ["id\teffLength", "g1\t100", "g2\t110"])

        logs = (out / "GeneCount_Assigned_logs.xls").read_text().splitlines()
        self.assertEqual(logs[0], "metric\tS1")
        self.assertIn("Assigned\t100", logs[1])

    def test_uneven_gene_sets_counts_filled_fpkm_warns(self):
        """S2 lacks g2: counts cell becomes 0, FPKM/TPM keep NaN + warn."""
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.count").write_text(count_tsv(["g1", "g2"]))
        (indir / "S2.count").write_text(count_tsv(["g1"]))

        proc = self._run("-i", str(indir), "-o", str(self.tmp / "out"))
        self.assertEqual(proc.returncode, 0, proc.stderr)

        counts = (self.tmp / "out" / "count.matrix.tsv").read_text().splitlines()
        self.assertEqual(counts[2], "g2\t11\t0")  # filled, integer
        self.assertIn("filled with 0", proc.stderr)

        fpkm = (self.tmp / "out" / "GeneExpression_FPKM.xls").read_text().splitlines()
        self.assertEqual(fpkm[2].split("\t")[1:], ["1.6", ""])  # NaN kept in S2
        self.assertIn("fpkm matrix keeps 1 missing values", proc.stderr)

    def test_efflen_mismatch_warns_first_file_wins(self):
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.count").write_text(count_tsv(["g1", "g2"]))
        (indir / "S2.count").write_text(
            count_tsv(["g1", "g2"], eff={"g2": 999}))  # annotation drift

        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("effLength differs for 1 shared gene(s)", proc.stderr)
        efflen = (indir / "effLength.txt").read_text().splitlines()
        self.assertEqual(efflen[2], "g2\t110")  # first file's value kept

    def test_native_fc_tsv_rejected_with_guidance(self):
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.fc.tsv").write_text(
            "# Program: featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tS1\n")
        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("run-featurecounts.R", proc.stderr)
        self.assertIn("--fc-script", proc.stderr)

    def test_fc_tsv_alongside_counts_is_ignored(self):
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.count").write_text(count_tsv(["g1"]))
        (indir / "S2.count").write_text(count_tsv(["g1"]))
        (indir / "stale.fc.tsv").write_text("Geneid\tChr\n")
        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Ignoring 1 native *.fc.tsv", proc.stderr)

    def test_duplicate_gene_ids_fail(self):
        indir = self.tmp / "in"
        indir.mkdir()
        body = count_tsv(["g1", "g1"])
        (indir / "S1.count").write_text(body)
        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Duplicate gene ids in", proc.stderr)

    def test_unreadable_file_reports_path(self):
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.count").write_bytes(b"\xff\xfe\x00binary")
        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("S1.count", proc.stderr)
        self.assertIn("Failed to read", proc.stderr)

    def test_missing_columns_fail(self):
        indir = self.tmp / "in"
        indir.mkdir()
        (indir / "S1.count").write_text("id\tcounts\ng1\t5\n")
        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Missing columns", proc.stderr)

    def test_empty_dir_fails(self):
        indir = self.tmp / "in"
        indir.mkdir()
        proc = self._run("-i", str(indir))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No *.count files", proc.stderr)

    def test_missing_input_dir_fails(self):
        proc = self._run("-i", str(self.tmp / "nope"))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not found", proc.stderr)


if __name__ == "__main__":
    unittest.main()
