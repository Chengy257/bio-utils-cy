"""Tests for bin/sra_prefetch_batch.sh: dry-run paths + stubbed e2e via real GNU parallel.

sra-tools are pinned via config/env.local.sh on this machine; tests that need
"tools missing" or stub tools override the BUC_* slots / PATH explicitly.
"""

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402

NO_TOOLS_ENV = {
    "PATH": "/usr/bin:/bin",
    "BUC_SKIP_LOCAL_ENV": "1",
    "BUC_PREFETCH_BIN": "",
    "BUC_FASTERQ_DUMP_BIN": "",
    "BUC_FASTQ_DUMP_BIN": "",
    "BUC_PARALLEL_BIN": "",
    "BUC_FASTQC_BIN": "",
}

PREFETCH_STUB = """#!/bin/bash
set -e
srr="$1"; cache="$3"
mkdir -p "${cache}/${srr}"
head -c 64 /dev/zero > "${cache}/${srr}/${srr}.sra"
"""

DUMP_STUB = """#!/bin/bash
set -e
sra="$1"; outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir) outdir="$2"; shift 2 ;;
        *) shift ;;
    esac
done
srr="$(basename "${sra}" .sra)"
echo "FASTQDATA" > "${outdir}/${srr}_1.fastq.gz"
"""


class TestSraPrefetchBatch(ScriptTestCase):

    def _run(self, *args):
        return self.run_script("sra_prefetch_batch.sh", *args)

    def _write_accs(self, text="SRR1\n# comment\n\nSRR2\r\nSRR1\n"):
        return self.write("accs.txt", text)

    # ------------------------------------------------------------------
    # Dry run
    # ------------------------------------------------------------------

    def test_dry_run_works_without_tools_installed(self):
        out = self.tmp / "out"
        accs = self._write_accs()
        with mock.patch.dict(os.environ, NO_TOOLS_ENV):
            proc = self._run("-i", accs, "-d", out, "--dry-run")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # plan on stdout, logs on stderr
        self.assertIn("[prefetch + dump] SRR1", proc.stdout)
        self.assertIn("[prefetch + dump] SRR2", proc.stdout)
        self.assertIn("DRY RUN plan", proc.stderr)
        self.assertIn("reference it by name", proc.stderr)  # soft resolver warned
        # nothing created
        self.assertFalse(out.exists())

    def test_dry_run_classifies_cache_and_done(self):
        out = self.tmp / "out"
        (out / "sra_cache" / "SRR3").mkdir(parents=True)
        (out / "sra_cache" / "SRR3" / "SRR3.sra").write_bytes(b"x" * 32)
        (out / "SRR2.fastq.gz").write_text("old")
        accs = self.write("accs.txt", "SRR1\nSRR2\nSRR3\n")
        with mock.patch.dict(os.environ, NO_TOOLS_ENV):
            proc = self._run("-i", accs, "-d", out, "--dry-run")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("convert cached .sra] SRR3", proc.stdout)
        self.assertIn("[skip] SRR2", proc.stderr)
        self.assertNotIn("SRR2", proc.stdout)

    def test_dry_run_warns_on_non_run_accessions(self):
        accs = self._write_accs("SRR1\nBADACC\n")
        with mock.patch.dict(os.environ, NO_TOOLS_ENV):
            proc = self._run("-i", accs, "-d", self.tmp / "out", "--dry-run")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("does not look like a run accession", proc.stderr)

    def test_invalid_max_size_rejected(self):
        accs = self._write_accs("SRR1\n")
        proc = self._run("-i", accs, "--max-size", "100X", "--dry-run")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Invalid --max-size", proc.stderr)

    # ------------------------------------------------------------------
    # End-to-end with stub prefetch/dump and real GNU parallel
    # ------------------------------------------------------------------

    def _make_stubs(self):
        stubdir = self.tmp / "stubbin"
        stubdir.mkdir()
        prefetch = stubdir / "prefetch"
        prefetch.write_text(PREFETCH_STUB)
        prefetch.chmod(0o755)
        dump = stubdir / "fasterq-dump"
        dump.write_text(DUMP_STUB)
        dump.chmod(0o755)
        return stubdir

    def _stub_env(self, stubdir):
        return {
            "PATH": f"{stubdir}:/usr/bin:/bin",
            "BUC_SKIP_LOCAL_ENV": "1",
            "BUC_PREFETCH_BIN": "",
            "BUC_FASTERQ_DUMP_BIN": "",
            "BUC_FASTQ_DUMP_BIN": "",
            "BUC_PARALLEL_BIN": "",
            "BUC_FASTQC_BIN": "",
        }

    def test_end_to_end_downloads_and_converts(self):
        stubdir = self._make_stubs()
        out = self.tmp / "dl"
        accs = self.write("accs.txt", "SRR1\nSRR2\n")
        env = self._stub_env(stubdir)
        env["BUC_FASTQ_DUMP_BIN"] = ""  # bulk mode uses fastq-dump; force stub dump
        # use --sc so the dump binary is fasterq-dump (our stub name), bulk flags fine
        with mock.patch.dict(os.environ, env):
            proc = self._run("-i", accs, "-d", out, "--sc", "-t", "2")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((out / "SRR1_1.fastq.gz").exists())
        self.assertTrue((out / "SRR2_1.fastq.gz").exists())
        self.assertEqual((out / "worker_logs" / "SRR1.rc").read_text(), "0")
        self.assertIn("ok: 2, failed: 0", proc.stderr)
        # cache populated by stub prefetch
        self.assertTrue((out / "sra_cache" / "SRR1" / "SRR1.sra").exists())
        # no staging leftovers
        self.assertEqual(list((out / ".staging").glob("*")), [])

    def test_end_to_exit_code_2_on_dump_failure(self):
        stubdir = self.tmp / "stubbin2"
        stubdir.mkdir()
        (stubdir / "prefetch").write_text(PREFETCH_STUB)
        (stubdir / "prefetch").chmod(0o755)
        # dump stub fails only for SRR2
        (stubdir / "fasterq-dump").write_text(
            DUMP_STUB.replace('echo "FASTQDATA"',
                              '[[ "$(basename "${sra}" .sra)" == SRR2 ]] && exit 1; echo "FASTQDATA"'))
        (stubdir / "fasterq-dump").chmod(0o755)

        out = self.tmp / "dl2"
        accs = self.write("accs.txt", "SRR1\nSRR2\n")
        with mock.patch.dict(os.environ, self._stub_env(stubdir)):
            proc = self._run("-i", accs, "-d", out, "--sc", "-t", "2")
        self.assertEqual(proc.returncode, 2, proc.stderr)
        self.assertTrue((out / "SRR1_1.fastq.gz").exists())
        self.assertFalse((out / "SRR2_1.fastq.gz").exists())
        self.assertEqual((out / "failed_accessions.txt").read_text(), "SRR2\n")
        self.assertIn("ok: 1, failed: 1", proc.stderr)

    def test_rerun_skips_finished_accessions(self):
        stubdir = self._make_stubs()
        out = self.tmp / "dl3"
        accs = self.write("accs.txt", "SRR1\n")
        with mock.patch.dict(os.environ, self._stub_env(stubdir)):
            first = self._run("-i", accs, "-d", out, "--sc")
            self.assertEqual(first.returncode, 0, first.stderr)
            # mark the output file with new content to prove it is not re-made
            (out / "SRR1_1.fastq.gz").write_text("PRESERVED")
            second = self._run("-i", accs, "-d", out, "--sc")
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("Already done: 1", second.stderr)
        self.assertEqual((out / "SRR1_1.fastq.gz").read_text(), "PRESERVED")


if __name__ == "__main__":
    unittest.main()
