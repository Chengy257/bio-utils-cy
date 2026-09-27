"""Offline tests for bin/search_sra_riboseq.py.

All NCBI Entrez traffic is mocked; no test performs network I/O. The
BioProject docsum fixture mirrors the real efetch structure captured from
NCBI (RecordSet > DocumentSummary > Project, with a container-style
<ProjectID> element -- the structure that made BioProjectID come out empty
before v1.1.0).
"""

import contextlib
import io
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

import search_sra_riboseq as srr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase


def load_fixture():
    return (Path(__file__).parent / "fixtures" / "bioproject_docsum.xml").read_bytes()


class _FakeHandle:
    """Duck-typed Entrez handle: Entrez.read() gets ._payload, .read() gets ._raw."""

    def __init__(self, payload=None, raw=None):
        self._payload = payload
        self._raw = raw

    def read(self):
        if self._raw is None:
            raise AssertionError("raw handle.read() not expected here")
        return self._raw

    def close(self):
        pass


class EntrezMock:
    """Patch Bio.Entrez esearch/efetch/esummary/read with canned responses.

    esearch dispatches on (db, term); efetch always returns the fixture XML;
    esummary returns the canned runs list.
    """

    def __init__(self, bioproject_ids, sra_ids_by_term, runs_items,
                 bioproject_xml=None):
        self.bioproject_ids = bioproject_ids
        self.sra_ids_by_term = sra_ids_by_term
        self.runs_items = runs_items
        self.bioproject_xml = load_fixture() if bioproject_xml is None else bioproject_xml
        self.esearch_calls = []  # (db, term)

    # -- fake Entrez functions ------------------------------------------
    def _esearch(self, db=None, term=None, **kw):
        self.esearch_calls.append((db, term))
        if db == "bioproject":
            return _FakeHandle(payload={"IdList": list(self.bioproject_ids)})
        if db == "sra":
            ids = self.sra_ids_by_term.get(term, [])
            return _FakeHandle(payload={"IdList": list(ids)})
        raise AssertionError(f"unexpected esearch db={db} term={term}")

    def _efetch(self, **kw):
        return _FakeHandle(raw=self.bioproject_xml)

    def _esummary(self, **kw):
        return _FakeHandle(payload=list(self.runs_items))

    @staticmethod
    def _read(handle):
        return handle._payload

    # -- context manager -------------------------------------------------
    def __enter__(self):
        p = mock.patch.object
        self._cms = [
            p(srr.Entrez, "esearch", side_effect=self._esearch),
            p(srr.Entrez, "efetch", side_effect=self._efetch),
            p(srr.Entrez, "esummary", side_effect=self._esummary),
            p(srr.Entrez, "read", side_effect=self._read),
        ]
        for cm in self._cms:
            cm.__enter__()
        return self

    def __exit__(self, *exc):
        for cm in reversed(self._cms):
            cm.__exit__(*exc)
        return False


RUNS_SHARED = [
    {"Runs": '<Run acc="SRR1" total_spots="100"/><Run acc="ERR9001"/>'},
    {"Runs": '<Run acc="DRR500"/>'},
    {"Runs": '<Run acc="SRR1"/>'},  # duplicate across items
]


class TestParseKeywords(ScriptTestCase):

    def test_classifies_by_content(self):
        result = srr.parse_keywords("ribosome profiling,translatome,polysome profiling")
        self.assertEqual(result["strict"], ["ribosome profiling"])
        self.assertEqual(result["medium"], ["translatome", "polysome profiling"])

    def test_case_insensitive_and_empty_tokens(self):
        result = srr.parse_keywords(" RPF-seq ,, TE profiling ,")
        self.assertEqual(result["strict"], ["RPF-seq"])
        self.assertEqual(result["medium"], ["TE profiling"])


class TestBuildQueries(ScriptTestCase):

    def test_taxid_adds_organism_filter(self):
        s = srr.NCBISearcher(email="t@example.com", api_delay=0)
        queries = s.build_queries(3702)
        self.assertEqual([name for name, _ in queries], ["STRICT", "MEDIUM"])
        for _, q in queries:
            self.assertIn(" AND txid3702[Organism]", q)
        strict_q = dict(queries)["STRICT"]
        self.assertIn('"ribosome profiling"[All Fields]', strict_q)

    def test_no_taxid_no_filter(self):
        s = srr.NCBISearcher(email="t@example.com", api_delay=0)
        for _, q in s.build_queries(None):
            self.assertNotIn("txid", q)


class TestParseBioprojectXml(ScriptTestCase):

    def test_container_projectid_falls_back_to_accession(self):
        root = srr.ET.fromstring(load_fixture())
        projects = root.findall(".//Project")
        self.assertEqual(len(projects), 2)
        rec = srr.NCBISearcher._parse_bioproject_xml(projects[0])
        self.assertEqual(rec["Accession"], "PRJNA111")
        # <ProjectID> is a container (text is whitespace) -> accession fallback
        self.assertEqual(rec["BioProjectID"], "PRJNA111")
        self.assertEqual(rec["Title"], "Ribo-seq study one")
        self.assertEqual(rec["Organism"], "Species one")
        self.assertIn("Ribosome profiling", rec["Description"])

    def test_projectlinks_outside_project_not_parsed_as_record(self):
        root = srr.ET.fromstring(load_fixture())
        # .//Project must match exactly the two record projects, not the
        # <ProjectIDRef>-bearing siblings
        self.assertEqual(len(root.findall(".//Project")), 2)

    def test_direct_text_projectid_preserved(self):
        xml = """<RecordSet><DocumentSummary uid="9">
        <Project><ProjectID>99999</ProjectID>
        <ProjectDescr><Title>T</Title></ProjectDescr></Project>
        </DocumentSummary></RecordSet>"""
        project = srr.ET.fromstring(xml).find(".//Project")
        rec = srr.NCBISearcher._parse_bioproject_xml(project)
        self.assertEqual(rec["BioProjectID"], "99999")
        self.assertEqual(rec["Accession"], "")


class TestLinkSra(ScriptTestCase):

    def test_extracts_err_and_drr_and_dedupes(self):
        s = srr.NCBISearcher(email="t@example.com", api_delay=0)
        with EntrezMock(["1", "2"], {"PRJNA111[BioProject]": ["1", "2"]}, RUNS_SHARED):
            runs = s.link_sra("PRJNA111")
        # ERR + DRR + SRR extracted; SRR1 duplicated in item 3 is deduped
        self.assertEqual(runs, ["DRR500", "ERR9001", "SRR1"])

    def test_no_hits_returns_empty(self):
        s = srr.NCBISearcher(email="t@example.com", api_delay=0)
        with EntrezMock([], {}, []):
            self.assertEqual(s.link_sra("PRJNA999"), [])


class TestSearchPipeline(ScriptTestCase):

    def _make_searcher(self, query_type="both"):
        return srr.NCBISearcher(
            email="t@example.com", api_delay=0, query_type=query_type)

    def _mock(self):
        return EntrezMock(
            bioproject_ids=["4000001", "4000002"],
            sra_ids_by_term={
                "PRJNA111[BioProject]": ["1", "2"],
                "PRJNA222[BioProject]": ["1"],  # overlapping run -> dedup
            },
            runs_items=RUNS_SHARED,
        )

    def test_end_to_end_mocked(self):
        out = self.tmp / "results"
        with self._mock() as m:
            result = self._make_searcher().search("3702", out)

        # dedup across projects: both projects report the same 3 runs
        self.assertEqual(result["srr_ids"], ["DRR500", "ERR9001", "SRR1"])

        # SRA queries used the ArchiveID accession, not numeric UIDs
        sra_terms = [term for db, term in m.esearch_calls if db == "sra"]
        self.assertEqual(sorted(set(sra_terms)),
                         ["PRJNA111[BioProject]", "PRJNA222[BioProject]"])
        self.assertNotIn("4000001[BioProject]", sra_terms)

        # output files
        srr_file = out / "srr_ids.txt"
        self.assertEqual(srr_file.read_text().split(),
                         ["DRR500", "ERR9001", "SRR1"])
        self.assertIn("report", result["output_files"])
        self.assertIn("Total SRR IDs extracted:     3",
                      (out / "search_report.txt").read_text())

        # TSV: one row per project, SRR_Count = 3, BioProjectID = accession
        tsv = (out / "bioprojects_summary.tsv").read_text().splitlines()
        self.assertEqual(tsv[0].split("\t")[0], "BioProjectID")
        rows = [line.split("\t") for line in tsv[1:]]
        self.assertEqual(len(rows), 2)
        self.assertEqual({r[0] for r in rows}, {"PRJNA111", "PRJNA222"})
        for r in rows:
            self.assertEqual(r[5], "3")

    def test_strict_tier_runs_single_bioproject_query(self):
        out = self.tmp / "results_strict"
        with self._mock() as m:
            self._make_searcher(query_type="strict").search("3702", out)
        bp_terms = [term for db, term in m.esearch_calls if db == "bioproject"]
        self.assertEqual(len(bp_terms), 1)

    def test_unresolvable_species_raises(self):
        s = self._make_searcher()
        with mock.patch.object(s, "species_to_taxid", return_value=None):
            with self.assertRaises(ValueError):
                s.search("Unknownus fictionalus", self.tmp)

    def test_empty_tier_raises(self):
        s = srr.NCBISearcher(
            email="t@example.com", api_delay=0, query_type="medium",
            custom_keywords={"strict": ["ribosome profiling"], "medium": []},
        )
        with self.assertRaises(ValueError):
            s.search("3702", self.tmp)

    def test_report_and_summary_share_one_timestamp(self):
        # statistics text is generated once and reused for file + stdout
        stats = srr.NCBISearcher._generate_statistics(
            [{"BioProjectID": "PRJX", "SRR_Count": 1, "Title": "t"}],
            ["SRR1"], "test")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            srr.NCBISearcher._print_summary(stats, {"report": "x"})
        printed = buf.getvalue()
        self.assertIn("Total BioProjects found:", printed)
        self.assertIn("Total SRR IDs extracted:", printed)


if __name__ == "__main__":
    unittest.main()
