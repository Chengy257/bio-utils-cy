"""Functional tests for bin/genome_format_converter.sh.

The kent tools are not installed on this machine, so end-to-end tests
run against stub executables that replicate the kent CLI argument
conventions. This verifies pipeline wiring, lazy per-mode tool
resolution, TMPDIR plumbing, and the staged (atomic) output without
needing the real tools.
"""

import os
import shutil
import subprocess
import unittest
from pathlib import Path

from common import REPO_ROOT, ScriptTestCase

HAS_KENT = all(
    shutil.which(t) for t in ("gtfToGenePred", "genePredToBed", "gff3ToGenePred", "genePredToGtf")
)

REAL_GTF = (
    'chr1\ttest\tgene\t100\t900\t.\t+\t.\tgene_id "g1";\n'
    'chr1\ttest\ttranscript\t100\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\ttest\texon\t100\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";\n'
    'chr1\ttest\texon\t500\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "2";\n'
    'chr1\ttest\tCDS\t100\t300\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\ttest\tCDS\t500\t700\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
)

STUB_GTF2GP = (
    '#!/bin/bash\n'
    '# gtfToGenePred INPUT stdout   (pipe stage: write to stdout)\n'
    '# gtfToGenePred INPUT OUTPUT   (gtf2gp: write to OUTPUT)\n'
    'if [[ "$2" == "stdout" ]]; then cat "$1"; else cat "$1" > "$2"; fi\n'
)
STUB_GFF3TOGP = '#!/bin/bash\n# gff3ToGenePred INPUT stdout\ncat "$1"\n'
STUB_GP2BED = '#!/bin/bash\n# genePredToBed stdin OUTPUT\ncat > "$2"\n'
STUB_GP2GTF = '#!/bin/bash\n# genePredToGtf file stdin OUTPUT\ncat > "$3"\n'
# Failing variant: writes output then exits non-zero (partial output scenario).
STUB_GP2BED_FAIL = '#!/bin/bash\ncat > "$2"\nexit 1\n'
# TMPDIR-aware variant for gtf2gp: verifies $TMPDIR is exported.
STUB_GTF2GP_TMP = (
    '#!/bin/bash\n'
    '[[ -n "${TMPDIR:-}" ]] || exit 9\n'
    'printf "STAGED" > "$2"\n'
)


class GenomeFormatConverterTest(ScriptTestCase):

    def setUp(self):
        super().setUp()
        self.stubs = self.tmp / "stubs"
        self.stubs.mkdir()
        for name, body in [
            ("gtfToGenePred", STUB_GTF2GP),
            ("gff3ToGenePred", STUB_GFF3TOGP),
            ("genePredToBed", STUB_GP2BED),
            ("genePredToGtf", STUB_GP2GTF),
        ]:
            p = self.stubs / name
            p.write_text(body)
            p.chmod(0o755)

    def _run(self, *args):
        return self.run_script("genome_format_converter.sh", *args)

    def _stub_args(self):
        return (
            "--gtftogenepred", self.stubs / "gtfToGenePred",
            "--genepredtobed", self.stubs / "genePredToBed",
            "--gff3togenepred", self.stubs / "gff3ToGenePred",
            "--genepredtogtf", self.stubs / "genePredToGtf",
        )

    def test_gtf2bed_pipeline_wiring(self):
        inp = self.write("in.gtf", "GPDATA_MARKER\n")
        out = self.tmp / "out.bed"
        proc = self._run("-i", inp, "-m", "gtf2bed", "-o", out, *self._stub_args())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "GPDATA_MARKER\n")

    def test_gff2bed_and_gff2gtf_wiring(self):
        inp = self.write("in.gff3", "GPDATAMARKER\n")
        out = self.tmp / "out.bed"
        proc = self._run("-i", inp, "-m", "gff2bed", "-o", out, *self._stub_args())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "GPDATAMARKER\n")
        out2 = self.tmp / "out.gtf"
        proc = self._run("-i", inp, "-m", "gff2gtf", "-o", out2, *self._stub_args())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out2.read_text(), "GPDATAMARKER\n")

    def test_auto_output_name_keeps_dotted_directory_and_last_ext(self):
        # Regression: ${INPUT%%.*} truncated at the first dot and lost
        # the directory; now only the last extension is replaced.
        sub = self.tmp / "v1.2"
        sub.mkdir()
        inp = sub / "anno.v2.gff3"
        inp.write_text("GPDATAMARKER\n")
        proc = self._run("-i", inp, "-m", "gff2bed", *self._stub_args())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((sub / "anno.v2.bed").exists(), proc.stderr)

    def test_lazy_tool_resolution(self):
        # gtf2bed must not require gff3ToGenePred/genePredToGtf.
        inp = self.write("in.gtf", "GPDATA\n")
        out = self.tmp / "out.bed"
        proc = self._run(
            "-i", inp, "-m", "gtf2bed", "-o", out,
            "--gtftogenepred", self.stubs / "gtfToGenePred",
            "--genepredtobed", self.stubs / "genePredToBed",
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_failed_pipeline_leaves_no_partial_output(self):
        inp = self.write("in.gtf", "GPDATA\n")
        out = self.tmp / "out.bed"
        fail_stub = self.stubs / "genePredToBedFail"
        fail_stub.write_text(STUB_GP2BED_FAIL)
        fail_stub.chmod(0o755)
        proc = self._run(
            "-i", inp, "-m", "gtf2bed", "-o", out,
            "--gtftogenepred", self.stubs / "gtfToGenePred",
            "--genepredtobed", fail_stub,
        )
        self.assertEqual(proc.returncode, 1)
        self.assertFalse(out.exists())  # staged output never moved into place
        leftovers = list(self.tmp.glob("gfc.*")) + list(Path("/tmp").glob("gfc.*"))
        self.assertEqual(leftovers, [])  # trap cleaned the staging file

    def test_tmpdir_option_is_used(self):
        inp = self.write("in.gtf", "x\n")
        out = self.tmp / "out.gp"
        tmpdir = self.tmp / "mytmp"
        tmpdir.mkdir()
        stub = self.stubs / "gtfToGenePredTmp"
        stub.write_text(STUB_GTF2GP_TMP)
        stub.chmod(0o755)
        proc = self._run("-i", inp, "-m", "gtf2gp", "-o", out, "-t", tmpdir,
                         "--gtftogenepred", stub)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "STAGED")

    def test_config_slot_resolution(self):
        # BUC_GTFTOGENEPRED_BIN (derived slot) must be honored without
        # an explicit CLI option.
        inp = self.write("in.gtf", "GPDATA\n")
        out = self.tmp / "out.gp"
        env = dict(os.environ)
        env["BUC_GTFTOGENEPRED_BIN"] = str(self.stubs / "gtfToGenePred")
        proc = subprocess.run(
            ["bash", str(REPO_ROOT / "bin" / "genome_format_converter.sh"),
             "-i", str(inp), "-m", "gtf2gp", "-o", str(out)],
            capture_output=True, text=True, env=env, timeout=120,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "GPDATA\n")

    def test_explicit_tool_not_executable_errors(self):
        inp = self.write("in.gtf", "x\n")
        proc = self._run("-i", inp, "-m", "gtf2gp",
                         "--gtftogenepred", "/nonexistent/tool")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("not executable", proc.stderr)

    def test_validation_errors(self):
        proc = self._run("-i", "/nonexistent.gtf", "-m", "gtf2bed")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Input file not found", proc.stderr)

        inp = self.write("in.gtf", "x\n")
        proc = self._run("-i", inp, "-m", "bogus")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Invalid mode", proc.stderr)

        proc = self._run("-i", inp, "-m", "gtf2bed", "-t", "/nonexistent/dir")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Temp directory does not exist", proc.stderr)

        proc = self._run("-i", inp, "-m", "gtf2bed", "--wat")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Unknown option", proc.stderr)

    def test_logs_go_to_stderr(self):
        inp = self.write("in.gtf", "GPDATA\n")
        out = self.tmp / "out.bed"
        proc = self._run("-i", inp, "-m", "gtf2bed", "-o", out, *self._stub_args())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "")


@unittest.skipUnless(HAS_KENT, "UCSC kent tools not on PATH")
class GenomeFormatConverterKentTest(ScriptTestCase):
    """Real conversions against the kent tools found on PATH."""

    def test_real_gtf2gp(self):
        inp = self.write("real.gtf", REAL_GTF)
        out = self.tmp / "real.gp"
        proc = self.run_script("genome_format_converter.sh", "-i", inp, "-m", "gtf2gp", "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        lines = out.read_text().splitlines()
        self.assertTrue(lines)
        fields = lines[0].split("\t")
        self.assertEqual(fields[0], "t1")       # name
        self.assertEqual(fields[1], "chr1")
        self.assertEqual(fields[2], "+")
        self.assertEqual(fields[3], "99")       # txStart, 0-based
        self.assertGreaterEqual(len(fields), 10)

    def test_real_gtf2bed(self):
        inp = self.write("real.gtf", REAL_GTF)
        out = self.tmp / "real.bed"
        proc = self.run_script("genome_format_converter.sh", "-i", inp, "-m", "gtf2bed", "-o", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        fields = out.read_text().splitlines()[0].split("\t")
        self.assertEqual(fields[0], "chr1")
        self.assertEqual(fields[3], "t1")
        self.assertEqual(fields[9], "2")        # blockCount: two exons
        # 1-based exon 100-300/500-900 -> 0-based BED
        self.assertEqual(fields[1], "99")


if __name__ == "__main__":
    unittest.main()
