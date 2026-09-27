"""Offline tests for bin/fetch_sra_metadata_ncbi.py.

All HTTP traffic is mocked via urllib.request.urlopen; no network is touched.
"""

import csv
import io
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fetch_sra_metadata_ncbi as fsn  # noqa: E402
from common import ScriptTestCase  # noqa: E402

RUNINFO_CSV = (
    "Run,ReleaseDate,LoadDate,spots,bases,BioSample,Experiment,LibraryName\n"
    "SRR1,2020-01-01,2020-01-02,100,200,SAMN1,SRX1,lib1\n"
    "SRR2,2020-01-01,2020-01-02,150,300,SAMN2,SRX2,lib2\n"
)

BIOSAMPLE_XML = """<?xml version="1.0" encoding="UTF-8"?>
<BioSampleSet>
  <BioSample accession="SAMN1" id="1">
    <Description><Title>sample one</Title>
      <Organism><OrganismName>Vibrio cholerae</OrganismName></Organism>
    </Description>
    <Attributes>
      <Attribute harmonized_name="strain" attribute_name="strain">O1</Attribute>
      <Attribute harmonized_name="isolate" attribute_name="isolate">rep1</Attribute>
    </Attributes>
  </BioSample>
  <BioSample accession="SAMN2" id="2">
    <Description><Title>sample two</Title></Description>
    <Attributes>
      <Attribute harmonized_name="strain" attribute_name="strain">O2</Attribute>
    </Attributes>
  </BioSample>
</BioSampleSet>
"""

EXPERIMENT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT><TITLE>exp one</TITLE></EXPERIMENT>
    <RUN_SET><RUN accession="SRR1" total_spots="100"/></RUN_SET>
    <LIBRARY_DESCRIPTOR>
      <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
      <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
      <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
      <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
    </LIBRARY_DESCRIPTOR>
    <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina HiSeq 4000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
  </EXPERIMENT_PACKAGE>
  <EXPERIMENT_PACKAGE>
    <EXPERIMENT><TITLE>exp two</TITLE></EXPERIMENT>
    <RUN_SET><RUN accession="SRR2" total_spots="150"/></RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""


def _fake_response(text):
    """A urlopen-shaped context manager with http.client-style headers."""
    resp = mock.Mock()
    resp.headers.get_content_charset.return_value = "utf-8"
    resp.read.return_value = text.encode()
    resp.__enter__ = mock.Mock(return_value=resp)
    resp.__exit__ = mock.Mock(return_value=False)
    return resp


def _make_urlopen(call_texts):
    """Return (fake_urlopen, urls_seen); call_texts consumed in order."""
    urls_seen = []

    def fake_urlopen(req, timeout=None):
        urls_seen.append(req.full_url)
        text = call_texts[min(len(urls_seen) - 1, len(call_texts) - 1)]
        return _fake_response(text)

    return fake_urlopen, urls_seen


class TestSanitizeUrl(ScriptTestCase):

    def test_api_key_masked(self):
        url = "https://x/efetch.fcgi?db=sra&api_key=SECRET123&id=SRR1"
        self.assertNotIn("SECRET123", fsn._sanitize_url(url))
        self.assertIn("api_key=***", fsn._sanitize_url(url))

    def test_url_without_key_unchanged(self):
        self.assertEqual(fsn._sanitize_url("https://x/?db=sra"),
                         "https://x/?db=sra")


class TestFetchUrlLogging(ScriptTestCase):

    def test_debug_log_masks_api_key(self):
        url = "https://eutils/efetch.fcgi?db=sra&api_key=TOPSECRET&id=1"
        with self.assertLogs(fsn.logger, level="DEBUG") as captured:
            with mock.patch.object(
                    fsn, "urlopen",
                    side_effect=AssertionError("no network in tests")):
                # urlopen raises, but the debug line must already be emitted
                with self.assertRaises((AssertionError, RuntimeError)):
                    fsn._fetch_url(url)
        joined = "\n".join(captured.output)
        self.assertNotIn("TOPSECRET", joined)
        self.assertIn("api_key=***", joined)


class TestRetryPolicy(ScriptTestCase):

    def test_4xx_not_retried(self):
        err = fsn.HTTPError("http://x", 404, "Not Found", None, None)
        with mock.patch.object(fsn, "urlopen", side_effect=err):
            with mock.patch.object(fsn.time, "sleep") as sleeper:
                with self.assertRaises(fsn._ClientHTTPError):
                    fsn._fetch_with_retry("http://x", max_retries=3)
        sleeper.assert_not_called()

    def test_429_retried(self):
        err = fsn.HTTPError("http://x", 429, "Too Many Requests", None, None)
        ok = mock.Mock()
        ok.headers.get_content_charset.return_value = "utf-8"
        ok.read.return_value = b"data"
        ok.__enter__ = mock.Mock(return_value=ok)
        ok.__exit__ = mock.Mock(return_value=False)
        responses = [err, ok]
        with mock.patch.object(fsn, "urlopen", side_effect=responses):
            with mock.patch.object(fsn.time, "sleep"):
                self.assertEqual(fsn._fetch_with_retry("http://x", max_retries=3),
                                 "data")

    def test_5xx_retried_then_exhausted(self):
        err = fsn.HTTPError("http://x", 503, "Service Unavailable", None, None)
        with mock.patch.object(fsn, "urlopen", side_effect=err):
            with mock.patch.object(fsn.time, "sleep"):
                with self.assertRaises(RuntimeError) as ctx:
                    fsn._fetch_with_retry("http://x", max_retries=2)
        self.assertIn("All 2 retries exhausted", str(ctx.exception))


class TestReadIds(ScriptTestCase):

    def test_dedupe_and_comments(self):
        f = self.write("ids.txt", "SRR1\n# c\n\nSRR1\nERR2\n")
        self.assertEqual(fsn._read_ids(str(f)), ["SRR1", "ERR2"])

    def test_warns_on_non_run_accessions(self):
        f = self.write("ids2.txt", "SRR1\nGARBAGE\n")
        with self.assertLogs(fsn.logger, level="WARNING"):
            fsn._read_ids(str(f))


class TestBiosampleKeyConsistency(ScriptTestCase):

    def test_extract_and_merge_share_keys(self):
        # a lowercase-keyed RunInfo row (older header style) must be both
        # extracted for fetching AND merged afterwards
        rows = [{"Run": "SRR9", "bio_sample": "SAMEX9"}]
        accs = fsn._extract_biosample_accessions(rows)
        self.assertEqual(accs, ["SAMEX9"])
        merged = fsn.merge_metadata(rows, {"SAMEX9": {"organism": "X"}}, {})
        self.assertEqual(merged[0].get("biosample_organism"), "X")

    def test_uppercase_key_still_merged(self):
        rows = [{"Run": "SRR1", "BioSample": "SAMN1"}]
        merged = fsn.merge_metadata(rows, {"SAMN1": {"strain": "O1"}}, {})
        self.assertEqual(merged[0]["biosample_strain"], "O1")


class TestParsers(ScriptTestCase):

    def test_fetch_biosample_xml_parses_attributes(self):
        fake, urls = _make_urlopen([BIOSAMPLE_XML])
        with mock.patch.object(fsn, "urlopen", side_effect=fake):
            data = fsn.fetch_biosample_xml(["SAMN1", "SAMN2"], email="a@b.c")
        self.assertEqual(data["SAMN1"]["strain"], "O1")
        self.assertEqual(data["SAMN1"]["organism"], "Vibrio cholerae")
        self.assertEqual(data["SAMN1"]["description_title"], "sample one")
        self.assertEqual(data["SAMN2"]["strain"], "O2")
        # tool + email sent with the request
        self.assertIn("email=a%40b.c", urls[0])
        self.assertIn("tool=fetch_sra_metadata_ncbi", urls[0])

    def test_fetch_experiment_xml_parses_fields(self):
        fake, urls = _make_urlopen([EXPERIMENT_XML])
        with mock.patch.object(fsn, "urlopen", side_effect=fake):
            data = fsn.fetch_experiment_xml(["SRR1", "SRR2"], email="a@b.c")
        self.assertEqual(data["SRR1"]["experiment_title"], "exp one")
        self.assertEqual(data["SRR1"]["library_strategy"], "RNA-Seq")
        self.assertEqual(data["SRR1"]["instrument_model"], "Illumina HiSeq 4000")
        self.assertEqual(data["SRR2"]["experiment_title"], "exp two")

    def test_fetch_runinfo_parses_csv(self):
        fake, urls = _make_urlopen([RUNINFO_CSV])
        with mock.patch.object(fsn, "urlopen", side_effect=fake):
            rows = fsn.fetch_runinfo(["SRR1", "SRR2"], email="a@b.c")
        self.assertEqual(rows[0]["Run"], "SRR1")
        self.assertEqual(rows[1]["BioSample"], "SAMN2")

    def test_api_key_param_sent_but_url_logged_masked(self):
        fake, urls = _make_urlopen([RUNINFO_CSV])
        with mock.patch.object(fsn, "urlopen", side_effect=fake):
            with self.assertLogs(fsn.logger, level="DEBUG") as captured:
                fsn.fetch_runinfo(["SRR1"], email="a@b.c", api_key="SECRET")
        self.assertIn("api_key=SECRET", urls[0])
        self.assertFalse(any("SECRET" in line for line in captured.output))


class TestMergeAndWrite(ScriptTestCase):

    def test_merge_full_pipeline(self):
        rows = [{"Run": "SRR1", "BioSample": "SAMN1", "spots": "100"}]
        biosample = {"SAMN1": {"organism": "Vibrio cholerae"}}
        experiment = {"SRR1": {"experiment_title": "exp one"}}
        merged = fsn.merge_metadata(rows, biosample, experiment)
        self.assertEqual(merged[0]["biosample_organism"], "Vibrio cholerae")
        self.assertEqual(merged[0]["experiment_experiment_title"], "exp one")

    def test_write_tsv(self):
        out = self.tmp / "out.tsv"
        fsn.write_tsv([{"Run": "SRR1", "extra": "e"}], str(out))
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "Run\textra")
        self.assertEqual(lines[1], "SRR1\te")

    def test_collect_fieldnames_first_seen_order(self):
        rows = [{"a": "1", "b": "2"}, {"c": "3", "a": "4"}]
        self.assertEqual(fsn.collect_fieldnames(rows), ["a", "b", "c"])


if __name__ == "__main__":
    unittest.main()
