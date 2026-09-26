"""Functional tests for bin/calculate_nucleotide_diversity.py.

plink/vcftools are pinned via config/env.local.sh (BUC_PLINK_BIN /
BUC_VCFTOOLS_BIN), so these tests build a tiny bfile from a synthetic
VCF and run real end-to-end conversions — including the BED-start
off-by-one regression.
"""

import os
import subprocess
import unittest
from pathlib import Path

from common import REPO_ROOT, ScriptTestCase


def _resolve(var, tool):
    cmd = f'. "{REPO_ROOT}/config/env.sh" 2>/dev/null; buc_resolve_bin {var} {tool}'
    proc = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 else None


PLINK = _resolve("BUC_PLINK_BIN", "plink")
VCFTOOLS = _resolve("BUC_VCFTOOLS_BIN", "vcftools")
HAS_TOOLS = bool(PLINK and VCFTOOLS)

# 1-based positions: 95, 96, 300, 400. POS 95/96 straddle a BED start
# boundary used in the off-by-one regression below.
VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t95\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t1/1\t1/1
chr1\t96\t.\tC\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t0/1
chr1\t300\t.\tG\tA\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0
chr1\t400\t.\tT\tC\t.\t.\t.\tGT\t0/1\t1/1\t1/1\t0/0
"""


def read_tsv(path):
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@unittest.skipUnless(HAS_TOOLS, "plink/vcftools not resolvable via config/env.sh")
class NucleotideDiversityTest(ScriptTestCase):

    @classmethod
    def setUpClass(cls):
        import tempfile
        cls._class_tmp = tempfile.TemporaryDirectory()
        cls.bprefix = str(Path(cls._class_tmp.name) / "bfile")
        vcf = Path(cls._class_tmp.name) / "in.vcf"
        vcf.write_text(VCF_TEXT)
        subprocess.run([PLINK, "--vcf", str(vcf), "--make-bed", "--out", cls.bprefix],
                       capture_output=True, check=True)

    @classmethod
    def tearDownClass(cls):
        cls._class_tmp.cleanup()

    def setUp(self):
        super().setUp()
        # Run the script against the class-wide bfile.
        self.bprefix = type(self).bprefix

    def _run(self, bed_text, extra=()):
        bed = self.write("regions.bed", bed_text)
        out = self.tmp / "pi.tsv"
        proc = self.run_script(
            "calculate_nucleotide_diversity.py", "-b", bed,
            "-p", self.bprefix, "-o", out,
            "--plink-bin", PLINK, "--vcftools-bin", VCFTOOLS,
            *extra,
        )
        return proc, out

    def test_basic_pi_and_zero_variant_interval(self):
        bed = "chr1\t94\t300\tregA\nchr1\t500\t600\treg_empty\n"
        proc, out = self._run(bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = {r["name"]: r for r in read_tsv(out)}
        # [94,300) covers 1-based 95..300: sites 95, 96 and 300 (the
        # monomorphic site 300 is kept by plink and dilutes the mean).
        self.assertEqual(int(rows["regA"]["snp_count"]), 3)
        pi_a = float(rows["regA"]["pi_value"])
        self.assertGreater(pi_a, 0.0)
        self.assertLess(pi_a, 1.0)
        # No variant sites in [500,600): real zero, not a failure.
        self.assertEqual(int(rows["reg_empty"]["snp_count"]), 0)
        self.assertEqual(float(rows["reg_empty"]["pi_value"]), 0.0)

    def test_bed_start_offbyone_regression(self):
        # BED [95,100) covers 1-based 96..100; the SNP at POS 95 must be
        # EXCLUDED (the old code passed the 0-based start verbatim and
        # included it).
        bed = "chr1\t95\t100\tboundary\n"
        proc, out = self._run(bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(int(rows[0]["snp_count"]), 1)  # POS 96 only
        # And BED [94,100) (1-based 95..100) must include both SNPs.
        proc, out = self._run("chr1\t94\t100\tboundary2\n")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(int(rows[0]["snp_count"]), 2)

    def test_bed6_columns_accepted(self):
        # Real BED6 (6 columns) crashed the old 4-name read_csv.
        bed = "chr1\t94\t300\tregA\t0\t+\n"
        proc, out = self._run(bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(rows[0]["name"], "regA")
        self.assertEqual(int(rows[0]["snp_count"]), 3)

    def test_bed3_columns_generated_name(self):
        proc, out = self._run("chr1\t94\t300\n")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(rows[0]["name"], "chr1:94-300")

    def test_duplicate_names_no_collision(self):
        bed = "chr1\t94\t300\tdup\nchr1\t390\t410\tdup\n"
        proc, out = self._run(bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(len(rows), 2)
        # Interval 1 has 3 sites (one monomorphic, pi-diluted); interval
        # 2 has the single mixed site pos400. Both 0<pi<1 and distinct.
        self.assertEqual(int(rows[0]["snp_count"]), 3)
        self.assertEqual(int(rows[1]["snp_count"]), 1)
        pis = [float(r["pi_value"]) for r in rows]
        self.assertTrue(all(0.0 < p < 1.0 for p in pis))
        self.assertNotAlmostEqual(pis[0], pis[1], places=3)

    def test_bad_vcftools_binary_yields_na_rows_and_exit_1(self):
        bed = "chr1\t94\t300\tregA\n"
        proc, out = self._run(bed, ["--vcftools-bin", "/bin/false"])
        self.assertEqual(proc.returncode, 1, proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(rows[0]["name"], "regA")
        self.assertEqual(rows[0]["snp_count"], "")   # NA, not 0
        self.assertEqual(rows[0]["pi_value"], "")    # NA, not 0.0
        self.assertIn("failed", proc.stderr)

    def test_unusable_tool_path_clean_error(self):
        bed = self.write("regions.bed", "chr1\t94\t300\tregA\n")
        proc = self.run_script("calculate_nucleotide_diversity.py", "-b", bed,
                               "-p", self.bprefix, "-o", self.tmp / "o.tsv",
                               "--plink-bin", "/nonexistent/plink")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not executable", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_malformed_bed_lines_skipped(self):
        bed = (
            "chr1\tXX\t300\tbad\n"
            "chr1\t300\t200\tbackwards\n"
            "chr1\t94\n"
            "# comment\n"
            "chr1\t94\t300\tok\n"
        )
        proc, out = self._run(bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("malformed", proc.stderr)
        rows = read_tsv(out)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["name"], "ok")

    def test_empty_bed_exits_1(self):
        proc, _ = self._run("# only comment\n")
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
