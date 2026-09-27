"""Offline tests for bin/pride_metadata.py.

The `parse` subcommand is fully offline (local JSON fixtures); the `fetch`
subcommand's HTTP layer is mocked.
"""

import contextlib
import io
import json
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pride_metadata as pm  # noqa: E402
from common import ScriptTestCase  # noqa: E402

PROJECT_JSON = {
    "accession": "PXD000001",
    "title": "Test project",
    "projectDescription": "A synthetic PRIDE project.",
    "publicationDate": "2020-01-01",
    "submissionDate": "2019-12-01",
    "submitters": [
        {"name": "Alice", "email": "a@x.io", "affiliation": "Uni X"},
        {"name": "Bob", "email": "b@x.io", "affiliation": "Uni Y"},
    ],
    "labPIs": [{"name": "Prof. X"}],
    "organisms": [{"name": "Homo sapiens"}, {"name": "Mus musculus"}],
    "tissues": [{"name": "liver"}],
    "instrumentNames": ["Orbitrap Fusion"],
    "quantificationMethods": [{"name": "Label free"}],
    "ptmNames": ["Phosphorylation"],
    "experimentTypes": [{"name": "Protein identification"}],
    "keywords": ["proteomics", "synthetic"],
    "doi": "10.9999/fake",
    "pubmedIds": ["12345678"],
    "numFiles": 3,
}


class TestExtractProjectMetadata(ScriptTestCase):

    def test_full_schema(self):
        row = pm.extract_project_metadata(PROJECT_JSON)
        self.assertEqual(row["accession"], "PXD000001")
        self.assertEqual(row["submitter_name"], "Alice")  # first submitter only
        self.assertEqual(row["lab_pi"], "Prof. X")
        self.assertEqual(row["species"], "Homo sapiens; Mus musculus")
        self.assertEqual(row["tissue"], "liver")
        self.assertEqual(row["instrument"], "Orbitrap Fusion")
        self.assertEqual(row["quantification_method"], "Label free")
        self.assertEqual(row["modification"], "Phosphorylation")
        self.assertEqual(row["keywords"], "proteomics; synthetic")
        self.assertEqual(row["pubmed_ids"], "12345678")
        self.assertEqual(row["num_files"], "3")

    def test_missing_fields_default_empty(self):
        row = pm.extract_project_metadata({"accession": "PXD1"})
        self.assertEqual(row["title"], "")
        self.assertEqual(row["species"], "")
        self.assertEqual(row["submitter_email"], "")


class TestExtractInstruments(ScriptTestCase):

    def test_instruments_cv_list_fallback(self):
        data = {"instruments": [{"name": "Q Exactive"}]}
        self.assertEqual(pm._extract_instruments(data), "Q Exactive")

    def test_additional_attributes_fallback(self):
        data = {"additionalAttributes": [
            {"cvParam": {"name": "instrument model", "value": "timsTOF"}}]}
        self.assertEqual(pm._extract_instruments(data), "timsTOF")

    def test_none_of_them(self):
        self.assertEqual(pm._extract_instruments({}), "")


class TestJoinCvList(ScriptTestCase):

    def test_dicts_plain_strings_and_garbage(self):
        self.assertEqual(pm._join_cv_list([{"name": "A"}, {"value": "B"}, "C"]), "A; B; C")
        self.assertEqual(pm._join_cv_list("not-a-list"), "")
        self.assertEqual(pm._join_cv_list([]), "")


class TestParseJsonFile(ScriptTestCase):

    def test_root_accession(self):
        f = self.write("p.json", json.dumps(PROJECT_JSON))
        row, err = pm.parse_json_file(Path(f))
        self.assertIsNone(err)
        self.assertEqual(row["accession"], "PXD000001")
        self.assertEqual(row["source_file"], "p.json")

    def test_wrapper_dict(self):
        f = self.write("w.json", json.dumps({"project": PROJECT_JSON}))
        row, err = pm.parse_json_file(Path(f))
        self.assertIsNone(err)
        self.assertEqual(row["accession"], "PXD000001")

    def test_projects_list_merged_with_separator(self):
        other = dict(PROJECT_JSON, accession="PXD000002", title="Second")
        f = self.write("m.json", json.dumps({"projects": [PROJECT_JSON, other]}))
        row, err = pm.parse_json_file(Path(f))
        self.assertIsNone(err)
        self.assertIn("PXD000001; PXD000002", row["accession"])
        self.assertIn("Test project; Second", row["title"])

    def test_invalid_json_returns_error(self):
        f = self.write("bad.json", "{not json")
        row, err = pm.parse_json_file(Path(f))
        self.assertIsNone(row)
        self.assertIn("Invalid JSON", err)

    def test_non_dict_returns_error(self):
        f = self.write("arr.json", json.dumps([1, 2]))
        row, err = pm.parse_json_file(Path(f))
        self.assertIsNone(row)
        self.assertIn("Expected dict", err)

    def test_unrecognized_dict_falls_back_on_title(self):
        f = self.write("t.json", json.dumps({"title": "Only a title"}))
        row, err = pm.parse_json_file(Path(f))
        self.assertIsNone(err)
        self.assertEqual(row["title"], "Only a title")


class TestFindJsonFiles(ScriptTestCase):

    def test_non_recursive_excludes_subdirs(self):
        good = self.write("a.json", "{}")
        self.write("b.txt", "nope")
        sub = self.tmp / "sub"
        sub.mkdir()
        (sub / "c.json").write_text("{}")
        found = pm.find_json_files([str(self.tmp)], recursive=False)
        self.assertEqual([p.name for p in found], ["a.json"])

    def test_recursive_includes_subdirs(self):
        self.write("a.json", "{}")
        sub = self.tmp / "sub"
        sub.mkdir()
        (sub / "c.json").write_text("{}")
        found = pm.find_json_files([str(self.tmp)], recursive=True)
        self.assertEqual(len(found), 2)

    def test_missing_path_warned_and_skipped(self):
        good = self.write("a.json", "{}")
        with self.assertLogs(pm.logger, level="WARNING"):
            found = pm.find_json_files([str(good), str(self.tmp / "nope")])
        self.assertEqual([p.name for p in found], ["a.json"])


class TestReadAccessions(ScriptTestCase):

    def test_dedupe_comments_and_warn(self):
        f = self.write("pxd.txt", "PXD1\n# c\n\nPXD1\nnotapxd\n")
        with self.assertLogs(pm.logger, level="WARNING"):
            accs = pm._read_accessions(str(f))
        self.assertEqual(accs, ["PXD1", "notapxd"])


class TestWriteTable(ScriptTestCase):

    def test_tsv_and_csv_files(self):
        rows = [{"a": "1", "b": "x"}]
        tsv = self.tmp / "o.tsv"
        pm._write_table(rows, str(tsv), fmt="tsv")
        self.assertEqual(tsv.read_text().splitlines(), ["a\tb", "1\tx"])
        csvf = self.tmp / "o.csv"
        pm._write_table(rows, str(csvf), fmt="csv")
        self.assertEqual(csvf.read_text().splitlines(), ["a,b", "1,x"])

    def test_stdout_csv_default_for_parse(self):
        rows = [{"a": "1"}]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            pm._write_table(rows, "-", fmt="csv")
        self.assertEqual(buf.getvalue().splitlines(), ["a", "1"])


class TestWriteXlsxFallback(ScriptTestCase):

    def test_falls_back_to_csv_without_openpyxl(self):
        rows = [{"a": "1", "b": "x"}]
        xlsx = self.tmp / "o.xlsx"
        with mock.patch.dict(sys.modules, {"openpyxl": None}):
            with self.assertLogs(pm.logger, level="WARNING"):
                pm.write_xlsx(rows, str(xlsx))
        csvf = self.tmp / "o.csv"
        self.assertTrue(csvf.exists())
        self.assertFalse(xlsx.exists())


class TestCmdParse(ScriptTestCase):

    def _run_parse(self, args):
        argv = ["pride_metadata.py", "parse"] + [str(a) for a in args]
        with mock.patch.object(sys, "argv", argv):
            try:
                pm.main()
                return 0
            except SystemExit as exc:
                return exc.code

    def test_mixed_files_partial_success(self):
        good = self.write("good.json", json.dumps(PROJECT_JSON))
        bad = self.write("bad.json", "{broken")
        out = self.tmp / "parsed.csv"
        with self.assertLogs(pm.logger, level="WARNING") as captured:
            code = self._run_parse([good, bad, "-o", out, "-f", "csv"])
        self.assertEqual(code, 0)
        body = out.read_text().splitlines()
        self.assertIn("PXD000001", body[1])
        self.assertTrue(any("had errors" in line for line in captured.output))

    def test_all_broken_exits_nonzero(self):
        bad = self.write("bad.json", "{broken")
        out = self.tmp / "parsed2.csv"
        code = self._run_parse([bad, "-o", out, "-f", "csv"])
        self.assertEqual(code, 1)
        self.assertFalse(out.exists())

    def test_no_json_files_found(self):
        code = self._run_parse([self.tmp / "empty_dir_xyz", "-o", "-"])
        self.assertEqual(code, 1)


class TestTopLevelVersion(ScriptTestCase):

    def test_version_prints(self):
        with mock.patch.object(sys, "argv", ["pride_metadata.py", "--version"]):
            with self.assertRaises(SystemExit) as ctx:
                pm.main()
        self.assertEqual(ctx.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
