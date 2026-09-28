"""Functional tests for bin/plot_qpcr_analysis.R.

Fixture layout mirrors a qPCR instrument export: 34 junk header rows,
one real column-header row, data rows (15 columns; sample=col4,
gene=col5, CT=col15) and 5 trailing rows.

Hand-computed expectations: ACTIN CT = 25 everywhere; GENE1 CT = 28 in
WT (control) replicates and 26 in TX (ddCT = (26-25) - (28-25) = -2 ->
relative expression 4.0); GENE2 CT = 30 in WT and TX except one
Undetermined well mapped to ct-max=45.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "ggplot2", "ggsci", "patchwork")


def build_qpcr_csv(samples, genes, ct_of, undetermined=()):
    """rows: (sample, gene, ct) repeated twice per pair unless in
    undetermined (then the second well reads 'Undetermined')."""
    rows = []
    for s in samples:
        for g in genes:
            for rep in range(2):
                ct = "Undetermined" if (s, g, rep) in undetermined else str(ct_of[(s, g)])
                row = [""] * 15
                row[3] = s
                row[4] = g
                row[14] = ct
                rows.append(",".join(row))
    lines = [",".join(f"hdr{i}_{j}" for j in range(15)) for i in range(34)]
    lines.append(",".join(f"col{j}" for j in range(15)))
    lines += rows
    lines += [",".join(f"tail{i}_{j}" for j in range(15)) for i in range(5)]
    return "\n".join(lines) + "\n"


def read_tsv(path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


@unittest.skipUnless(RS, "R with getopt+ggplot2+ggsci+patchwork not available")
class TestPlotQpcrAnalysis(ScriptTestCase):
    def setUp(self):
        super().setUp()
        ct_of = {}
        for s in ("WT_rep1", "WT_rep2"):
            ct_of[(s, "ACTIN")] = 25.0
            ct_of[(s, "GENE1")] = 28.0
            ct_of[(s, "GENE2")] = 30.0
        ct_of[("TX1_rep1", "ACTIN")] = 25.0
        ct_of[("TX1_rep1", "GENE1")] = 26.0
        ct_of[("TX1_rep1", "GENE2")] = 30.0
        self.csv = self.write("qpcr.csv", build_qpcr_csv(
            ["WT_rep1", "WT_rep2", "TX1_rep1"],
            ["ACTIN", "GENE1", "GENE2"],
            ct_of,
            undetermined={("TX1_rep1", "GENE2", 1)}))

    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # outputs are CWD-relative
        try:
            return subprocess.run([RS, str(BIN / "plot_qpcr_analysis.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def _table(self, prefix="qp"):
        return read_tsv(self.tmp / f"{prefix}_parsed_result_table.tsv")

    def test_ddct_values_and_control_normalisation(self):
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT", "-o", "qp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = self._table()
        g1 = [r for r in rows if r["gene"] == "GENE1"]
        # control group (both WT replicates) mean deltaCT = 3.0
        wt = [float(r["rel_exp"]) for r in g1 if r["sample"].startswith("WT")]
        for v in wt:
            self.assertAlmostEqual(v, 1.0, places=6)
        tx = [float(r["rel_exp"]) for r in g1 if r["sample"] == "TX1_rep1"]
        for v in tx:
            self.assertAlmostEqual(v, 4.0, places=6)  # 2^-(-2)
        # no row pairing: control rows carry the group mean
        for r in g1:
            self.assertAlmostEqual(float(r["control_mean_deltaCT"]), 3.0, places=6)

    def test_undetermined_mapped_to_ct_max(self):
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT", "-o", "qp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = [r for r in self._table()
                if r["gene"] == "GENE2" and r["sample"] == "TX1_rep1"]
        # 2 target wells x 2 reference wells -> combination rows; the
        # distinct target CTs are what matters
        cts = sorted({float(r["CT"]) for r in rows})
        self.assertEqual(cts, [30.0, 45.0])  # second well was Undetermined
        rels = sorted({round(float(r["rel_exp"]), 8) for r in rows})
        self.assertEqual(len(rels), 2)
        self.assertAlmostEqual(rels[0], 2.0 ** -15, places=8)
        self.assertAlmostEqual(rels[1], 1.0, places=6)

    def test_control_prefix_group_and_exact(self):
        # group prefix: "WT" -> both WT replicates form the control group
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT", "-o", "qp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Control samples: WT_rep1, WT_rep2", proc.stderr)
        # exact single sample also works
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT_rep1", "-o", "qp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Control samples: WT_rep1", proc.stderr)

    def test_control_partial_match_rejected(self):
        # "T" is a substring of every sample but matches neither exactly
        # nor as a prefix; the old grep silently accepted it
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "T", "-o", "qp")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Control sample 'T' not found", proc.stderr)

    def test_value_options_accepted(self):
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT", "-o", "qp",
                         "--skip-rows", "35", "--tail-rows", "5",
                         "--ct-max", "45", "--plot-width", "8",
                         "--plot-height", "8")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "qp_parsed_result_barplot.pdf").exists())
        self.assertGreater((self.tmp / "qp_parsed_result_barplot.pdf").stat().st_size,
                           3000)

    def test_bare_option_rejected(self):
        # regression: a bare flag was read as TRUE (read.csv(skip=TRUE))
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT",
                         "--skip-rows")
        self.assertEqual(proc.returncode, 1)

    def test_well_count_mismatch_warns_per_gene(self):
        text = build_qpcr_csv(["WT_rep1", "TX1_rep1"],
                              ["ACTIN", "GENE1"],
                              {("WT_rep1", "ACTIN"): 25.0,
                               ("WT_rep1", "GENE1"): 28.0,
                               ("TX1_rep1", "ACTIN"): 25.0,
                               ("TX1_rep1", "GENE1"): 26.0}).splitlines()
        # append one extra GENE1 well for TX1_rep1 (3 target vs 2 reference)
        extra = [""] * 15
        extra[3] = "TX1_rep1"
        extra[4] = "GENE1"
        extra[14] = "26.2"
        text.insert(34 + 1 + 8, ",".join(extra))  # after data, before tail rows
        self.write("qpcr3.csv", "\n".join(text) + "\n")

        proc = self._run("-d", "qpcr3.csv", "-r", "ACTIN", "-c", "WT_rep1",
                         "-o", "qp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("TX1_rep1/GENE1 (3 target vs 2 reference)", proc.stderr)

    def test_missing_reference_gene(self):
        proc = self._run("-d", self.csv, "-r", "GAPDH", "-c", "WT", "-o", "qp")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Reference gene 'GAPDH' not found", proc.stderr)

    def test_table_is_tsv_with_named_columns(self):
        proc = self._run("-d", self.csv, "-r", "ACTIN", "-c", "WT", "-o", "qp")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        tsv = self.tmp / "qp_parsed_result_table.tsv"
        self.assertTrue(tsv.exists())
        self.assertFalse((self.tmp / "qp_parsed_result_table.xls").exists())
        expected = ["gene", "sample", "CT", "ref_CT", "deltaCT",
                    "control_mean_deltaCT", "deltadeltaCT", "rel_exp"]
        self.assertEqual(tsv.read_text().splitlines()[0].split("\t"), expected)

    def test_missing_required_args(self):
        proc = self._run("-d", self.csv)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Missing required arguments", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage", proc.stdout)


if __name__ == "__main__":
    unittest.main()
