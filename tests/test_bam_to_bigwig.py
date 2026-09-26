"""Functional tests for bin/bam_to_bigwig.sh.

Validation paths run without external tools (input validation happens
before tool resolution). End-to-end conversion tests build a tiny
synthetic BAM with the pinned samtools and run the pinned bamCoverage;
they are skipped when either tool is not resolvable.
"""

import subprocess
import time
import unittest
from pathlib import Path

from common import REPO_ROOT, ScriptTestCase


def _resolve(var, tool):
    cmd = f'. "{REPO_ROOT}/config/env.sh" 2>/dev/null; buc_resolve_bin {var} {tool}'
    proc = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 else None


SAMTOOLS = _resolve("BUC_SAMTOOLS_BIN", "samtools")
BAMCOV = _resolve("BUC_BAMCOVERAGE_BIN", "bamCoverage")
HAS_TOOLS = bool(SAMTOOLS and BAMCOV)

SAM_TEXT = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "r1\t0\tchr1\t101\t60\t4M\t*\t0\t0\tACGT\t*\n"
    "r2\t0\tchr1\t201\t60\t4M\t*\t0\t0\tACGT\t*\n"
    "r3\t0\tchr1\t301\t60\t4M\t*\t0\t0\tACGT\t*\n"
)


def make_bam(tmpdir, name):
    """Create an indexed BAM from SAM_TEXT; returns the BAM path."""
    sam = Path(tmpdir) / f"{name}.sam"
    bam = Path(tmpdir) / f"{name}.bam"
    sam.write_text(SAM_TEXT)
    with open(bam, "wb") as out:
        subprocess.run([SAMTOOLS, "view", "-b", str(sam)], stdout=out, check=True)
    subprocess.run([SAMTOOLS, "index", str(bam)], check=True)
    return bam


class BamToBigwigValidationTest(ScriptTestCase):
    """Argument/input validation; no external tools required."""

    def test_missing_input(self):
        proc = self.run_script("bam_to_bigwig.sh")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Missing required option", proc.stderr)

    def test_invalid_norm(self):
        proc = self.run_script("bam_to_bigwig.sh", "-i", "x.bam", "-n", "BAD")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Invalid normalization", proc.stderr)

    def test_rpgc_requires_genome_size(self):
        proc = self.run_script("bam_to_bigwig.sh", "-i", "x.bam", "-n", "RPGC")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("requires --genome-size", proc.stderr)

    def test_rpgc_genome_size_must_be_integer(self):
        proc = self.run_script("bam_to_bigwig.sh", "-i", "x.bam", "-n", "RPGC", "-g", "abc")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("positive integer", proc.stderr)

    def test_numeric_options_validated(self):
        for flag in ("-t", "-p", "-b"):
            proc = self.run_script("bam_to_bigwig.sh", "-i", "x.bam", flag, "0")
            self.assertEqual(proc.returncode, 1, flag)
            self.assertIn("positive integer", proc.stderr)

    def test_list_with_comments_and_only_missing_bams(self):
        # The former first-line heuristic treated such a file as a single
        # BAM; now the list is parsed and the missing paths reported.
        listing = self.write("list.txt", "# comment\n\n/nonexistent/a.bam\n/nonexistent/b.bam\n")
        proc = self.run_script("bam_to_bigwig.sh", "-i", listing)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("BAM not found", proc.stderr)
        self.assertIn("No valid BAM files", proc.stderr)

    def test_nonexistent_input_errors(self):
        proc = self.run_script("bam_to_bigwig.sh", "-i", "/nonexistent/x.txt")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Input file not found", proc.stderr)


@unittest.skipUnless(HAS_TOOLS, "samtools/bamCoverage not resolvable via config/env.sh")
class BamToBigwigE2ETest(ScriptTestCase):
    """End-to-end conversions against the pinned deepTools/samtools."""

    def _run(self, *args):
        return self.run_script("bam_to_bigwig.sh", *args)

    def test_convert_single_bam_from_list_file(self):
        bam = make_bam(self.tmp, "s1")
        listing = self.write("list.txt", f"# sample list\n\n{bam}\n")
        proc = self._run("-i", listing, "-t", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        bw = self.tmp / "s1.bw"
        self.assertTrue(bw.exists())
        self.assertGreater(bw.stat().st_size, 100)
        self.assertIn("1/1 files converted", proc.stderr)

    def test_norm_none_is_explicit(self):
        # Regression: -n None used to omit the flag, so bamCoverage
        # silently applied its RPKM default; and on bash 4.2 the empty
        # argument array crashed under set -u.
        bam = make_bam(self.tmp, "s2")
        proc = self._run("-i", str(bam), "-n", "None", "-t", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "s2.bw").exists())

    def test_parallel_conversion(self):
        b1 = make_bam(self.tmp, "p1")
        b2 = make_bam(self.tmp, "p2")
        listing = self.write("list.txt", f"{b1}\n{b2}\n")
        outdir = self.tmp / "bw"
        proc = self._run("-i", listing, "-o", outdir, "-p", "2", "-t", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((outdir / "p1.bw").exists())
        self.assertTrue((outdir / "p2.bw").exists())
        self.assertIn("2/2 files converted", proc.stderr)

    def test_corrupt_bam_fails_batch_without_partial_output(self):
        good = make_bam(self.tmp, "good")
        bad = self.tmp / "bad.bam"
        bad.write_text("this is not a BAM\n")
        listing = self.write("list.txt", f"{good}\n{bad}\n")
        proc = self._run("-i", listing, "-t", "1")
        self.assertEqual(proc.returncode, 1, proc.stderr)
        self.assertIn("1/2 files converted", proc.stderr)
        self.assertFalse((self.tmp / "bad.bw").exists())  # no partial output

    def test_index_lifecycle(self):
        bam = make_bam(self.tmp, "stale")
        bai = Path(str(bam) + ".bai")

        bai.unlink()  # missing index -> must be created
        proc = self._run("-i", str(bam), "-n", "None", "-t", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Indexing:", proc.stderr)
        m_after_index = bai.stat().st_mtime

        time.sleep(1.1)
        bam.write_bytes(bam.read_bytes())  # BAM now newer than index (stale)
        proc = self._run("-i", str(bam), "-n", "None", "-f", "-t", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Indexing:", proc.stderr)  # stale index rebuilt
        self.assertGreater(bai.stat().st_mtime, m_after_index)

        time.sleep(1.1)
        proc = self._run("-i", str(bam), "-n", "None", "-f", "-t", "1")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertNotIn("Indexing:", proc.stderr)  # fresh index kept


if __name__ == "__main__":
    unittest.main()
