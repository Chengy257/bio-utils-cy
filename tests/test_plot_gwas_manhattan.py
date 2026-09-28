"""Functional tests for bin/plot_gwas_manhattan.R.

Manhattan + QQ paths run against real qqman with synthetic EMMAX/map
data. The regional path needs Gviz/magick (not installed) and is only
statically reviewed.
"""

import os
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase, find_r_pair  # noqa: E402

BIN = Path(__file__).resolve().parent.parent / "bin"

RS, RLIBS = find_r_pair("getopt", "ggplot2", "qqman", "data.table")


def build_gwas(n_per_chr=(("1", 200), ("2", 150)), extra_rows=None, seed=42):
    import random
    random.seed(seed)
    emmax, mapr, i = [], [], 0
    for c, n in n_per_chr:
        pos = sorted(random.sample(range(1, 10_000_000), n))
        for p in pos:
            i += 1
            snp = f"snp{i:05d}"
            pv = min(1.0, max(1e-12, random.random() ** 3))
            emmax.append(f"{snp}\t0.1\t0.05\t{pv:.3e}")
            mapr.append(f"{c}\t{snp}\t8\t{p}")
    for row in (extra_rows or []):
        emmax.append(row[0])
        mapr.append(row[1])
    return "\n".join(emmax) + "\n", "\n".join(mapr) + "\n"


@unittest.skipUnless(RS, "R with getopt+ggplot2+qqman+data.table not available")
class TestPlotGwasManhattan(ScriptTestCase):
    def _run(self, *args):
        env = dict(os.environ)
        if RLIBS:
            env["R_LIBS"] = RLIBS
        old_cwd = os.getcwd()
        os.chdir(self.tmp)  # outputs are CWD-relative
        try:
            return subprocess.run([RS, str(BIN / "plot_gwas_manhattan.R")]
                                  + [str(a) for a in args],
                                  capture_output=True, text=True, timeout=600,
                                  env=env)
        finally:
            os.chdir(old_cwd)

    def test_full_run_produces_both_pdfs(self):
        emmax, mapr = build_gwas()
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        proc = self._run("-e", e, "-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        for name in ("res_manhattan.pdf", "res_qq.pdf"):
            p = self.tmp / name
            self.assertTrue(p.exists() and p.stat().st_size > 3000, name)
        self.assertIn("Total SNPs: 350", proc.stderr)

    def test_nonnumeric_chromosomes_dropped_with_warning(self):
        # regression: chrX/scaffold SNPs crashed qqman ("'times' invalid")
        emmax, mapr = build_gwas(extra_rows=[
            ("xsnp1\t0.2\t0.02\t1e-9", "chrX\txsnp1\t8\t500000"),
            ("xsnp2\t0.2\t0.02\t1e-12", "scaffold1\txsnp2\t8\t501000"),
        ])
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        proc = self._run("-e", e, "-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("non-numeric chromosomes", proc.stderr)
        self.assertIn("scaffold1", proc.stderr)
        self.assertIn("Total SNPs: 350", proc.stderr)  # only numeric chrs
        self.assertTrue((self.tmp / "res_qq.pdf").exists())

    def test_out_of_range_pvalues_dropped(self):
        emmax, mapr = build_gwas(extra_rows=[
            ("pz\t0.1\t0.05\t0.0", "1\tpz\t8\t20000000"),
            ("pn\t0.1\t0.05\t-1e-3", "2\tpn\t8\t20000001"),
            ("pb\t0.1\t0.05\t1.5", "2\tpb\t8\t20000002"),
        ])
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        proc = self._run("-e", e, "-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("out-of-range p-values", proc.stderr)
        self.assertIn("Total SNPs: 350", proc.stderr)

    def test_merge_drops_reported(self):
        emmax, mapr = build_gwas()
        # remove the last 10 SNPs from the map
        mapr = "\n".join(mapr.splitlines()[:-10]) + "\n"
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        proc = self._run("-e", e, "-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("10 EMMAX SNP(s) not present in the map file", proc.stderr)
        self.assertIn("Total SNPs: 340", proc.stderr)

    def test_chr_prefixes_normalised(self):
        # "chr01" style prefixes must map onto chromosome 1
        emmax, mapr = build_gwas(n_per_chr=(("chr01", 120), ("Chr002", 80)))
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        proc = self._run("-e", e, "-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Total SNPs: 200", proc.stderr)
        self.assertNotIn("non-numeric", proc.stderr)

    def test_bad_thresholds_rejected(self):
        emmax, mapr = build_gwas()
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        for flag, val in (("--signif-line", "0"), ("--suggest-line", "2")):
            proc = self._run("-e", e, "-m", m, "-o", "res", flag, val)
            self.assertEqual(proc.returncode, 1, flag)
            self.assertIn("(0, 1]", proc.stderr)

    def test_missing_inputs(self):
        emmax, mapr = build_gwas()
        e = self.write("gw.emmax", emmax)
        m = self.write("gw.map", mapr)
        proc = self._run("-e", self.tmp / "nope", "-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("EMMAX file not found", proc.stderr)
        proc = self._run("-m", m, "-o", "res")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Missing required argument: -e", proc.stderr)

    def test_help(self):
        proc = self._run("-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)


if __name__ == "__main__":
    unittest.main()
