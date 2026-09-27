"""Offline tests for bin/fetch_sra_metadata_ena.py.

urlopen is mocked everywhere; no network is touched.
"""

import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fetch_sra_metadata_ena as fen  # noqa: E402
from common import ScriptTestCase  # noqa: E402

TSV_HEAD = "run_accession\texperiment_title\n"
FIELDS = "run_accession,experiment_title"


def _tsv(rows):
    return TSV_HEAD + "".join(f"{acc}\ttitle-{acc}\n" for acc in rows)


class _Resp:
    def __init__(self, text):
        self.headers = mock.Mock()
        self.headers.get_content_charset.return_value = "utf-8"
        self._text = text.encode()

    def read(self):
        return self._text

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestBuildUrl(ScriptTestCase):

    def test_url_contains_encoded_params(self):
        url = fen.build_url(["SRR1", "ERR2"], "run_accession,fastq_ftp")
        self.assertIn("accession=SRR1%2CERR2", url)
        self.assertIn("result=read_run", url)
        self.assertIn("format=tsv", url)


class TestChunkList(ScriptTestCase):

    def test_chunking(self):
        items = [str(i) for i in range(7)]
        chunks = fen.chunk_list(items, 3)
        self.assertEqual([len(c) for c in chunks], [3, 3, 1])


class TestReadAccessions(ScriptTestCase):

    def test_dedupe_comments_blank_crlf(self):
        f = self.write("a.txt", "SRR1\n# c\n\nSRR2\r\nSRR1\n")
        self.assertEqual(fen.read_accessions(str(f)), ["SRR1", "SRR2"])

    def test_warns_on_non_run_accessions(self):
        f = self.write("b.txt", "SRR1\nGARBAGE\n")
        with self.assertLogs(fen.logger, level="WARNING"):
            fen.read_accessions(str(f))


class TestFetchTsv(ScriptTestCase):

    def test_404_raises_dedicated_error_without_retry(self):
        err = fen.HTTPError("http://x", 404, "Not Found", None, None)
        urlopen = mock.Mock(side_effect=err)
        with mock.patch.object(fen, "urlopen", urlopen):
            with self.assertRaises(fen.EnaNotFoundError):
                fen.fetch_tsv("http://x", retries=3)
        self.assertEqual(urlopen.call_count, 1)  # deterministic: no retry

    def test_500_retries_then_raises(self):
        err = fen.HTTPError("http://x", 500, "Boom", None, None)
        with mock.patch.object(fen, "urlopen", side_effect=err):
            with mock.patch.object(fen.time, "sleep"):
                with self.assertRaises(RuntimeError):
                    fen.fetch_tsv("http://x", retries=2)

    def test_success_returns_text(self):
        with mock.patch.object(fen, "urlopen", return_value=_Resp("a\tb\n1\t2\n")):
            self.assertEqual(fen.fetch_tsv("http://x"), "a\tb\n1\t2\n")


class TestFetchChunkRows(ScriptTestCase):

    def _urlopen_by_accession(self, bad_accs):
        """Return rows for good accessions; 404 when a batch contains any bad one."""

        def fake_urlopen(req, timeout=None):
            url = req.full_url
            accs = url.split("accession=")[1].split("&")[0].split("%2C")
            bad = [a for a in accs if a in bad_accs]
            if bad:
                raise fen.HTTPError(url, 404, "Not Found", None, None)
            return _Resp(_tsv(accs))

        return fake_urlopen

    def test_404_bisects_to_isolate_bad_accession(self):
        failed = []
        urlopen = self._urlopen_by_accession({"ERRBAD"})
        with mock.patch.object(fen, "urlopen", side_effect=urlopen):
            rows = fen.fetch_chunk_rows(
                ["SRR1", "SRR2", "ERRBAD", "SRR4"], FIELDS, 60, 3, failed)
        self.assertEqual(failed, ["ERRBAD"])
        self.assertEqual([r["run_accession"] for r in rows],
                         ["SRR1", "SRR2", "SRR4"])

    def test_single_accession_404_recorded(self):
        failed = []
        urlopen = self._urlopen_by_accession({"ERRBAD"})
        with mock.patch.object(fen, "urlopen", side_effect=urlopen):
            rows = fen.fetch_chunk_rows(["ERRBAD"], FIELDS, 60, 3, failed)
        self.assertEqual(rows, [])
        self.assertEqual(failed, ["ERRBAD"])

    def test_non_tsv_error_body_recorded(self):
        failed = []

        def fake_urlopen(req, timeout=None):
            return _Resp("The following accessions are not embl: SRR1")

        with mock.patch.object(fen, "urlopen", side_effect=fake_urlopen):
            rows = fen.fetch_chunk_rows(["SRR1", "SRR2"], FIELDS, 60, 3, failed)
        self.assertEqual(rows, [])
        self.assertEqual(failed, ["SRR1", "SRR2"])

    def test_single_field_response_without_tabs_is_kept(self):
        failed = []
        with mock.patch.object(fen, "urlopen",
                               return_value=_Resp("run_accession\nSRR1\n")):
            rows = fen.fetch_chunk_rows(["SRR1"], "run_accession", 60, 3, failed)
        self.assertEqual(rows, [{"run_accession": "SRR1"}])
        self.assertEqual(failed, [])


class TestMain(ScriptTestCase):

    def _run_main(self, accs_text, bad_accs, output):
        accfile = self.write("accs.txt", accs_text)

        def fake_urlopen(req, timeout=None):
            url = req.full_url
            accs = url.split("accession=")[1].split("&")[0].split("%2C")
            good = [a for a in accs if a not in bad_accs]
            if len(good) < len(accs):  # ENA 404s the whole batch if any accession is unknown
                raise fen.HTTPError(url, 404, "Not Found", None, None)
            return _Resp(_tsv(good))

        argv = ["fetch_sra_metadata_ena.py", "-i", str(accfile), "-o", str(output)]
        with mock.patch.object(fen, "urlopen", side_effect=fake_urlopen):
            with mock.patch.object(fen.time, "sleep"):
                with mock.patch.object(sys, "argv", argv):
                    try:
                        fen.main()
                        code = 0
                    except SystemExit as exc:
                        code = exc.code
        return code

    def test_partial_failure_writes_failed_list_and_exits_nonzero(self):
        out = self.tmp / "meta.tsv"
        code = self._run_main("SRR1\nSRR2\nERRBAD\nSRR4\n", {"ERRBAD"}, out)
        self.assertEqual(code, 1)
        body = out.read_text().splitlines()
        self.assertEqual(body[0], "run_accession\texperiment_title")
        accs = [line.split("\t")[0] for line in body[1:]]
        self.assertEqual(sorted(accs), ["SRR1", "SRR2", "SRR4"])
        failed_file = Path(str(out) + ".failed.txt")
        self.assertEqual(failed_file.read_text(), "ERRBAD\n")

    def test_all_good_exits_zero_without_failed_file(self):
        out = self.tmp / "meta2.tsv"
        code = self._run_main("SRR1\nSRR2\n", set(), out)
        self.assertEqual(code, 0)
        self.assertFalse(Path(str(out) + ".failed.txt").exists())

    def test_fieldnames_are_union_across_batches(self):
        # batch 1 (SRR1) returns cols A+B; batch 2 (SRR2) returns A+C
        accfile = self.write("accs3.txt", "SRR1\nSRR2\n")
        out = self.tmp / "union.tsv"

        def fake_urlopen(req, timeout=None):
            url = req.full_url
            if "SRR1" in url:
                return _Resp("run_accession\tb\nSRR1\t1\n")
            return _Resp("run_accession\tc\nSRR2\t2\n")

        argv = ["fetch_sra_metadata_ena.py", "-i", str(accfile),
                "-o", str(out), "--fields", "run_accession,b,c",
                "--batch-size", "1"]  # force two separate batches
        with mock.patch.object(fen, "urlopen", side_effect=fake_urlopen):
            with mock.patch.object(fen.time, "sleep"):
                with mock.patch.object(sys, "argv", argv):
                    fen.main()  # success: returns without SystemExit
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0].split("\t"),
                         ["run_accession", "b", "c"])
        self.assertEqual(len(lines), 3)

    def test_all_failed_exits_nonzero(self):
        out = self.tmp / "meta3.tsv"
        code = self._run_main("ERRBAD1\nERRBAD2\n",
                              {"ERRBAD1", "ERRBAD2"}, out)
        self.assertEqual(code, 1)
        self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
