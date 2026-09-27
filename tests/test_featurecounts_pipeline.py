"""End-to-end tests for bin/featurecounts_pipeline.sh.

Strandedness detection is exercised against stub infer_experiment.py output
captured verbatim from the real RSeQC tool on this machine (PE: "1++,1--,
2+-,2-+" / "1+-,1-+,2++,2--"; SE: "++,--" / "+-,-+"). The -s mapping was
verified empirically against real featureCounts on a synthetic PE BAM:
read1 on the gene's sense strand is counted by -s 1, not -s 2.
"""

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402

PE_FWD = 'Fraction of reads explained by "1++,1--,2+-,2-+": {v}'
PE_REV = 'Fraction of reads explained by "1+-,1-+,2++,2--": {v}'
SE_FWD = 'Fraction of reads explained by "++,--": {v}'
SE_REV = 'Fraction of reads explained by "+-,-+": {v}'


class TestFeaturecountsPipeline(ScriptTestCase):

    def setUp(self):
        super().setUp()
        self.stubs = self.tmp / "stubs"
        self.stubs.mkdir()
        self.bam_dir = self.tmp / "bam"
        self.bam_dir.mkdir()
        self.bed = self.write("genes.bed", "chr1\t0\t200\tgeneA\t0\t+\n")
        self.gtf = self.write(
            "genes.gtf",
            'chr1\ttest\texon\t1\t200\t.\t+\t.\tgene_id "geneA";\n',
        )

        samtools = self.stubs / "samtools"
        samtools.write_text(
            "#!/bin/bash\n"
            'if [[ "$1" == index ]]; then : > "$4.bai"; exit 0; fi\n'
            'if [[ "$1" == view ]]; then echo "${STUB_PAIRED:-0}"; exit 0; fi\n'
            "exit 1\n"
        )
        samtools.chmod(0o755)

        infer = self.stubs / "infer_experiment.py"
        infer.write_text(
            "#!/bin/bash\n"
            'if [[ -n "${STUB_FAIL:-}" ]]; then echo "fatal: boom" >&2; exit 1; fi\n'
            'cat "${STUB_LINES_FILE:?}"\n'
        )
        infer.chmod(0o755)

        fc = self.stubs / "featureCounts"
        fc.write_text(
            "#!/bin/bash\n"
            'printf "%s\\n" "$@" > "${STUB_FC_ARGS:?}"\n'
            'out=""; prev=""\n'
            'for a in "$@"; do [[ "$prev" == "-o" ]] && out="$a"; prev="$a"; done\n'
            'if [[ -n "${STUB_EMPTY_OUT:-}" ]]; then : > "$out"; else\n'
            '    printf "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\tS1\\ngeneA\\tchr1\\t1\\t200\\t+\\t200\\t42\\n" > "$out"\n'
            "fi\n"
        )
        fc.chmod(0o755)

        rscript = self.stubs / "Rscript"
        rscript.write_text('#!/bin/bash\nprintf "%s\\n" "$@" > "${STUB_RS_ARGS:?}"\n')
        rscript.chmod(0o755)

        self.fc_args = self.tmp / "fc_args.txt"
        self.rs_args = self.tmp / "rs_args.txt"

    def _env(self, **extra):
        env = {
            "BUC_SKIP_LOCAL_ENV": "1",
            "BUC_SAMTOOLS_BIN": str(self.stubs / "samtools"),
            "BUC_INFER_EXP_BIN": str(self.stubs / "infer_experiment.py"),
            "BUC_FEATURECOUNTS_BIN": str(self.stubs / "featureCounts"),
            "BUC_RSCRIPT_BIN": str(self.stubs / "Rscript"),
            "STUB_FC_ARGS": str(self.fc_args),
            "STUB_RS_ARGS": str(self.rs_args),
        }
        env.update(extra)
        return env

    def _stub_lines(self, *lines):
        f = self.tmp / "stub_lines.txt"
        f.write_text("\n".join(lines) + "\n")
        return f

    def _add_bam(self, name="S1_Aligned.sortedByCoord.out.bam"):
        (self.bam_dir / name).write_bytes(b"")

    def _run(self, *args, **extra_env):
        with mock.patch.dict(os.environ, self._env(**extra_env)):
            return self.run_script("featurecounts_pipeline.sh", *args)

    def _base(self):
        return ("-b", str(self.bam_dir), "-g", str(self.gtf), "-r", str(self.bed),
                "-o", str(self.tmp / "out"), "-t", "4")

    def _strand_arg(self, proc):
        args = self.fc_args.read_text().splitlines()
        return args[args.index("-s") + 1]

    # --- strandedness detection -> featureCounts -s mapping ---

    def test_pe_forward_maps_to_s1(self):
        """read1 on the gene's sense strand (secondstrand) -> -s 1."""
        self._add_bam()
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "This is PairEnd Data",
            PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("secondstrand", proc.stderr)
        self.assertEqual(self._strand_arg(proc), "1")

    def test_pe_reverse_maps_to_s2(self):
        """read1 antisense (firststrand/dUTP) -> -s 2."""
        self._add_bam()
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "This is PairEnd Data",
            PE_FWD.format(v="0.0900"), PE_REV.format(v="0.9000"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("firststrand", proc.stderr)
        self.assertEqual(self._strand_arg(proc), "2")

    def test_pe_mixed_is_unstranded(self):
        """0.45/0.48 split: neither orientation leads by > 0.3 -> -s 0."""
        self._add_bam()
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "This is PairEnd Data",
            PE_FWD.format(v="0.4500"), PE_REV.format(v="0.4800"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("unstranded", proc.stderr)
        self.assertEqual(self._strand_arg(proc), "0")

    def test_se_forward_and_reverse(self):
        self._add_bam()
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "This is SingleEnd Data",
            SE_FWD.format(v="0.9200"), SE_REV.format(v="0.0600"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self._strand_arg(proc), "1")

        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "This is SingleEnd Data",
            SE_FWD.format(v="0.0500"), SE_REV.format(v="0.9300"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self._strand_arg(proc), "2")

    def test_unrecognized_output_defaults_unstranded(self):
        self._add_bam()
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "some unexpected output", "no fractions here")))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self._strand_arg(proc), "0")

    def test_infer_failure_reports_stderr_and_defaults(self):
        """Tool failure must surface the real stderr, then default to -s 0."""
        self._add_bam()
        proc = self._run(*self._base(), STUB_FAIL="1",
                         STUB_LINES_FILE=str(self._stub_lines("ignored")))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("defaulting to unstranded", proc.stderr)
        self.assertIn("fatal: boom", proc.stderr)
        self.assertEqual(self._strand_arg(proc), "0")

    # --- PE/SE detection ---

    def test_pe_flag_follows_paired_reads(self):
        self._add_bam()
        lines = str(self._stub_lines("This is PairEnd Data",
                                     PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500")))
        proc = self._run(*self._base(), STUB_PAIRED="500", STUB_LINES_FILE=lines)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("-p", self.fc_args.read_text().splitlines())
        self.assertIn("Paired-end", proc.stderr)

        proc = self._run(*self._base(), STUB_PAIRED="0", STUB_LINES_FILE=lines)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertNotIn("-p", self.fc_args.read_text().splitlines())

    # --- tool resolution: CLI > BUC_*_BIN > PATH ---

    def test_cli_samtools_respected(self):
        self._add_bam()
        other = self.stubs / "samtools2"
        other.write_text('#!/bin/bash\nprintf "%s\\n" "$@" > "${STUB_ALT_ARGS:?}"\n'
                         'if [[ "$1" == index ]]; then : > "$4.bai"; exit 0; fi\n'
                         'if [[ "$1" == view ]]; then echo 0; exit 0; fi\n'
                         "exit 1\n")
        other.chmod(0o755)
        alt = self.tmp / "alt_args.txt"
        proc = self._run(*self._base(), "--samtools", str(other),
                         STUB_LINES_FILE=str(self._stub_lines("This is PairEnd Data",
                                                              PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))),
                         STUB_ALT_ARGS=str(alt))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("view", alt.read_text())

    def test_cli_infer_exp_respected(self):
        """CLI --infer-exp wins over BUC_INFER_EXP_BIN: the stub's own
        (different) fractions must drive the -s value."""
        self._add_bam()
        other = self.stubs / "infer2.py"
        # stub2 reports a REVERSE-stranded library; the env-pinned stub1
        # reports forward — so -s 2 proves the CLI value was used.
        other.write_text('#!/bin/bash\ncat "${STUB_LINES_FILE:?}"\n')
        other.chmod(0o755)
        proc = self._run(*self._base(), "--infer-exp", str(other),
                         STUB_LINES_FILE=str(self._stub_lines(
                             "This is PairEnd Data",
                             PE_FWD.format(v="0.0500"), PE_REV.format(v="0.9400"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("firststrand", proc.stderr)
        self.assertEqual(self._strand_arg(proc), "2")

    # --- sample naming from --bam-pattern ---

    def test_bam_pattern_derives_sample_name(self):
        self._add_bam("S2_star.bam")
        proc = self._run(*self._base(), "--bam-pattern", "*_star.bam",
                         STUB_LINES_FILE=str(self._stub_lines(
                             "This is PairEnd Data",
                             PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Processing: S2", proc.stderr)
        self.assertTrue((self.tmp / "out" / "expression" / "S2.fc.tsv").exists())

    def test_bam_pattern_equals_form(self):
        self._add_bam("S3_star.bam")
        proc = self._run(*self._base(), "--bam-pattern=*_star.bam",
                         STUB_LINES_FILE=str(self._stub_lines(
                             "This is PairEnd Data",
                             PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "out" / "expression" / "S3.fc.tsv").exists())

    # --- --fc-script branch ---

    def test_fc_script_receives_strand_and_prefix(self):
        self._add_bam()
        script = self.write("run_fc.R", "# stub\n")
        proc = self._run(*self._base(), "--fc-script", str(script),
                         STUB_LINES_FILE=str(self._stub_lines(
                             "This is PairEnd Data",
                             PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        args = self.rs_args.read_text().splitlines()
        self.assertEqual(len(args), 6)
        self.assertEqual(args[0], str(script))
        self.assertTrue(args[1].endswith("S1_Aligned.sortedByCoord.out.bam"))
        self.assertTrue(args[2].endswith("genes.gtf"))       # gtf
        self.assertEqual(args[3], "1")                       # detected strand
        self.assertEqual(args[4], "4")                       # threads
        self.assertTrue(args[5].endswith("/expression/S1"))  # output prefix
        self.assertFalse(self.fc_args.exists())              # featureCounts not called

    def test_fc_script_missing_fails_early(self):
        self._add_bam()
        proc = self._run(*self._base(), "--fc-script", str(self.tmp / "nope.R"),
                         STUB_LINES_FILE=str(self._stub_lines("x")))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("--fc-script not found", proc.stderr)

    # --- post-run verification & error paths ---

    def test_empty_featurecounts_output_fails(self):
        self._add_bam()
        proc = self._run(*self._base(), STUB_EMPTY_OUT="1",
                         STUB_LINES_FILE=str(self._stub_lines(
                             "This is PairEnd Data",
                             PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("produced no output", proc.stderr)

    def test_no_bam_files_fails(self):
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines("x")))
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No BAM files", proc.stderr)

    def test_threads_must_be_numeric(self):
        self._add_bam()
        for bad in ("abc", "0", "-2"):
            proc = self._run("-b", str(self.bam_dir), "-g", str(self.gtf),
                             "-r", str(self.bed), "-o", str(self.tmp / "out"),
                             "-t", bad,
                             STUB_LINES_FILE=str(self._stub_lines("x")))
            self.assertEqual(proc.returncode, 1, bad)
            self.assertIn("Threads", proc.stderr)

    def test_logs_go_to_stderr(self):
        self._add_bam()
        proc = self._run(*self._base(), STUB_LINES_FILE=str(self._stub_lines(
            "This is PairEnd Data",
            PE_FWD.format(v="0.9500"), PE_REV.format(v="0.0500"))))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "")
        self.assertIn("[INFO]", proc.stderr)


if __name__ == "__main__":
    unittest.main()
