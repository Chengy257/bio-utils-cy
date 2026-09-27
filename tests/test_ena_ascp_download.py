"""Offline tests for bin/ena_ascp_download.py.

urllib and subprocess are mocked everywhere; no network or ascp is touched.
"""

import json
import os
import stat
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import ena_ascp_download as ead  # noqa: E402
from common import ScriptTestCase  # noqa: E402


class _FakeResp:
    def __init__(self, data):
        self._data = data

    def read(self):
        return json.dumps(self._data).encode()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestFindAscp(ScriptTestCase):

    def test_explicit_path(self):
        fake = self.tmp / "ascp"
        fake.write_text("#!/bin/sh\n")
        fake.chmod(0o755)
        self.assertEqual(ead.find_ascp(str(fake)), str(fake.resolve()))

    def test_explicit_non_executable_raises(self):
        fake = self.tmp / "ascp"
        fake.write_text("data")
        with self.assertRaises(FileNotFoundError):
            ead.find_ascp(str(fake))

    def test_env_slot_used(self):
        fake = self.tmp / "ascp"
        fake.write_text("#!/bin/sh\n")
        fake.chmod(0o755)
        with mock.patch.dict(os.environ, {"BUC_ASCP_BIN": str(fake)}):
            self.assertEqual(ead.find_ascp(None), str(fake.resolve()))

    def test_path_lookup(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("BUC_ASCP_BIN", None)
            with mock.patch.object(ead.shutil, "which", return_value="/fake/path/ascp"):
                self.assertEqual(ead.find_ascp(None), "/fake/path/ascp")


class TestFindAscpKey(ScriptTestCase):

    def test_finds_key_in_etc(self):
        bindir = self.tmp / "bin"
        bindir.mkdir()
        (bindir / "ascp").write_text("#!/bin/sh\n")
        etc = self.tmp / "etc"
        etc.mkdir()
        key = etc / "aspera_id_rsa"
        key.write_text("KEY")
        self.assertEqual(ead.find_ascp_key(str(bindir / "ascp")), str(key))

    def test_missing_key_raises(self):
        bindir = self.tmp / "bin"
        bindir.mkdir()
        (bindir / "ascp").write_text("#!/bin/sh\n")
        with self.assertRaises(FileNotFoundError):
            ead.find_ascp_key(str(bindir / "ascp"))


class TestReadAccessions(ScriptTestCase):

    def test_file_dedupe_comments_blanks(self):
        f = self.write("accs.txt", "ERR1\n# comment\n\nERR1\nERR2\n")
        self.assertEqual(ead.read_accessions(str(f)), ["ERR1", "ERR2"])

    def test_comma_string_dedupe(self):
        self.assertEqual(ead.read_accessions("SRR1, SRR2,SRR1"),
                         ["SRR1", "SRR2"])

    def test_warns_on_non_run_accessions(self):
        with self.assertLogs(ead.logger, level="WARNING"):
            ead.read_accessions("not_an_accession")


class TestParseAsperaUrls(ScriptTestCase):

    def test_split_and_strip(self):
        self.assertEqual(
            ead.parse_aspera_urls(" host:/a/b_1.fastq.gz;host:/a/b_2.fastq.gz;;"),
            ["host:/a/b_1.fastq.gz", "host:/a/b_2.fastq.gz"])


class TestResolveEnaUrls(ScriptTestCase):

    def _entries(self, n):
        return [
            {"run_accession": f"ERR{i:04d}",
             "fastq_aspera": f"fasp.sra.ebi.ac.uk:/vol1/ERR{i:04d}_1.fastq.gz",
             "fastq_bytes": "1234"}
            for i in range(n)
        ]

    def test_chunked_requests(self):
        urls_seen = []

        def fake_urlopen(req, timeout=None):
            urls_seen.append(req.full_url)
            return _FakeResp(self._entries(2))

        with mock.patch.object(ead.urllib.request, "urlopen", side_effect=fake_urlopen):
            records = ead.resolve_ena_urls([f"ERR{i:04d}" for i in range(1050)],
                                           retries=1)

        self.assertEqual(len(urls_seen), 3)  # 500 + 500 + 50
        counts = [u.split("accession=")[1].split("&")[0].count(",") + 1
                  for u in urls_seen]
        self.assertEqual(counts, [500, 500, 50])
        # ftp field is no longer requested; bytes field is
        self.assertNotIn("fastq_ftp", urls_seen[0])
        self.assertIn("fastq_bytes", urls_seen[0])
        self.assertEqual(len(records), 6)
        self.assertEqual(records[0], ("ERR0000",
                                      "fasp.sra.ebi.ac.uk:/vol1/ERR0000_1.fastq.gz",
                                      "1234"))

    def test_retry_then_success(self):
        responses = [ead.urllib.error.URLError("transient"), _FakeResp(self._entries(1))]
        with mock.patch.object(ead.urllib.request, "urlopen",
                               side_effect=responses):
            with mock.patch.object(ead.time, "sleep"):
                records = ead.resolve_ena_urls(["ERR0001"], retries=3)
        self.assertEqual(len(records), 1)

    def test_all_retries_exhausted_raises_runtimeerror(self):
        with mock.patch.object(ead.urllib.request, "urlopen",
                               side_effect=ead.urllib.error.URLError("down")):
            with mock.patch.object(ead.time, "sleep"):
                with self.assertRaises(RuntimeError):
                    ead.resolve_ena_urls(["ERR0001"], retries=2)

    def test_missing_aspera_skipped(self):
        entries = [{"run_accession": "ERR0001", "fastq_aspera": ""},
                   {"run_accession": "ERR0002", "fastq_aspera": "host:/x.gz"}]
        with mock.patch.object(ead.urllib.request, "urlopen",
                               return_value=_FakeResp(entries)):
            records = ead.resolve_ena_urls(["ERR0001", "ERR0002"], retries=1)
        self.assertEqual([r[0] for r in records], ["ERR0002"])


class TestRunAscp(ScriptTestCase):

    def _run(self, url, outdir, expected_size=None, timeout=3600, retries=1,
             subprocess_side_effect=None):
        captured = []

        def default_run(cmd, **kw):
            captured.append((cmd, kw))
            # simulate a successful transfer creating the file
            dest = os.path.join(outdir, url.rsplit("/", 1)[-1])
            Path(dest).write_bytes(b"x" * (expected_size or 3))
            return mock.Mock(returncode=0, stdout=b"", stderr=b"")

        se = subprocess_side_effect or default_run
        with mock.patch.object(ead.subprocess, "run", side_effect=se):
            dest = ead.run_ascp(
                url=url, outdir=str(outdir), ascp_bin="/usr/bin/ascp",
                ascp_key=str(self.write("key.pem", "K")), host="era-fasp@fasp.sra.ebi.ac.uk:",
                bandwidth="500M", timeout=timeout, retries=retries,
                expected_size=expected_size)
        return dest, captured

    def test_skips_existing_file_with_matching_size(self):
        outdir = self.tmp / "out"
        outdir.mkdir()
        (outdir / "f_1.fastq.gz").write_bytes(b"x" * 100)
        with mock.patch.object(ead.subprocess, "run",
                               side_effect=AssertionError("ascp should not run")):
            dest = ead.run_ascp(
                url="host:/v1/f_1.fastq.gz", outdir=str(outdir),
                ascp_bin="ascp", ascp_key="k", host="era-fasp@fasp.sra.ebi.ac.uk:",
                bandwidth="1M", expected_size=100)
        self.assertEqual(dest, str(outdir / "f_1.fastq.gz"))

    def test_truncated_existing_file_redownloaded(self):
        outdir = self.tmp / "out2"
        outdir.mkdir()
        (outdir / "f_1.fastq.gz").write_bytes(b"x" * 50)
        dest, captured = self._run("host:/v1/f_1.fastq.gz", outdir,
                                   expected_size=100)
        self.assertEqual(len(captured), 1)
        self.assertEqual(Path(dest).stat().st_size, 100)

    def test_existing_file_without_size_still_skipped(self):
        outdir = self.tmp / "out3"
        outdir.mkdir()
        (outdir / "f_1.fastq.gz").write_bytes(b"x" * 50)
        with mock.patch.object(ead.subprocess, "run",
                               side_effect=AssertionError("ascp should not run")):
            dest = ead.run_ascp(
                url="host:/v1/f_1.fastq.gz", outdir=str(outdir),
                ascp_bin="ascp", ascp_key="k", host="h:",
                bandwidth="1M", expected_size=None)
        self.assertTrue(dest)

    def test_host_prefix_applied_to_userless_url(self):
        outdir = self.tmp / "out4"
        outdir.mkdir()
        dest, captured = self._run("fasp.sra.ebi.ac.uk:/vol1/f.fastq.gz", outdir)
        cmd = captured[0][0]
        self.assertIn("era-fasp@fasp.sra.ebi.ac.uk:/vol1/f.fastq.gz", cmd)

    def test_url_with_user_used_verbatim(self):
        outdir = self.tmp / "out5"
        outdir.mkdir()
        dest, captured = self._run("otheruser@otherhost:/vol1/g.fastq.gz", outdir)
        cmd = captured[0][0]
        self.assertIn("otheruser@otherhost:/vol1/g.fastq.gz", cmd)

    def test_timeout_zero_disables_cap(self):
        outdir = self.tmp / "out6"
        outdir.mkdir()
        _, captured = self._run("host:/v1/f.gz", outdir, timeout=0)
        self.assertIsNone(captured[0][1]["timeout"])

    def test_retry_until_success(self):
        outdir = self.tmp / "out7"
        outdir.mkdir()
        calls = []

        def flaky(cmd, **kw):
            calls.append(cmd)
            if len(calls) < 3:
                raise ead.subprocess.TimeoutExpired(cmd, 10)
            Path(os.path.join(str(outdir), "f.gz")).write_bytes(b"xxx")
            return mock.Mock(returncode=0, stdout=b"", stderr=b"")

        with mock.patch.object(ead.time, "sleep"):
            dest, _ = self._run("host:/v1/f.gz", outdir, retries=3,
                                subprocess_side_effect=flaky)
        self.assertEqual(len(calls), 3)
        self.assertEqual(Path(dest).stat().st_size, 3)

    def test_all_retries_fail_raises(self):
        outdir = self.tmp / "out8"
        outdir.mkdir()
        err = ead.subprocess.CalledProcessError(1, "ascp", stderr=b"boom")

        def always_fail(cmd, **kw):
            raise err

        with mock.patch.object(ead.time, "sleep"):
            with self.assertRaises(RuntimeError) as ctx:
                self._run("host:/v1/f.gz", outdir, retries=2,
                          subprocess_side_effect=always_fail)
        self.assertIn("All 2 attempts failed", str(ctx.exception))


class TestMainExitCode(ScriptTestCase):

    def _setup_and_run(self, run_ascp_side_effects):
        accfile = self.write("accs.txt", "ERR0001\n")
        outdir = self.tmp / "dl"
        # fake ascp + key so tool resolution succeeds
        ascp = self.tmp / "ascp"
        ascp.write_text("#!/bin/sh\n")
        ascp.chmod(0o755)
        key = self.write("key.pem", "K")

        entries = [{"run_accession": "ERR0001",
                    "fastq_aspera": "fasp.sra.ebi.ac.uk:/vol1/a_1.fastq.gz;"
                                    "fasp.sra.ebi.ac.uk:/vol1/a_2.fastq.gz",
                    "fastq_bytes": "100;200"}]

        def fake_urlopen(req, timeout=None):
            return _FakeResp(entries)

        def fake_run_ascp(**kw):
            effect = run_ascp_side_effects.pop(0)
            filename = kw["url"].rsplit("/", 1)[-1]
            if isinstance(effect, Exception):
                raise effect
            dest = Path(str(outdir)) / filename
            dest.write_bytes(b"x" * (kw["expected_size"] or 1))
            return str(dest)

        argv = ["ena_ascp_download.py", "-i", str(accfile), "-o", str(outdir),
                "--ascp-bin", str(ascp), "--ascp-key", str(key)]
        with mock.patch.object(ead.urllib.request, "urlopen",
                               side_effect=fake_urlopen):
            with mock.patch.object(ead, "run_ascp", side_effect=fake_run_ascp):
                with mock.patch.object(sys, "argv", argv):
                    try:
                        ead.main()
                        exit_code = 0
                    except SystemExit as exc:
                        exit_code = exc.code
        return outdir, exit_code

    def test_partial_failure_exits_nonzero(self):
        outdir, code = self._setup_and_run([None, RuntimeError("net down")])
        self.assertEqual(code, 1)
        faillog = outdir / "failed_downloads.txt"
        self.assertTrue(faillog.exists())
        self.assertIn("a_2.fastq.gz", faillog.read_text())

    def test_all_success_exits_zero(self):
        outdir, code = self._setup_and_run([None, None])
        self.assertEqual(code, 0)
        self.assertFalse((outdir / "failed_downloads.txt").exists())


if __name__ == "__main__":
    unittest.main()
