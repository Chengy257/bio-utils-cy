"""Functional tests for bin/format_sra_summary.py."""

import csv

from common import ScriptTestCase


def write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    lines.extend("\t".join(r) for r in rows)
    path.write_text("\n".join(lines) + "\n")


class FormatSraSummaryTest(ScriptTestCase):
    def run_tool(self, inp, extra=()):
        out = self.tmp / "summary.tsv"
        proc = self.run_script("format_sra_summary.py", "--input", inp,
                               "--output", out, *extra)
        return proc, out

    def read_rows(self, out):
        with out.open() as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_field_remap_and_seqtype(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "LibraryStrategy", "Organism", "Instrument"],
                  [["SRR1", "WGS", "Oryza sativa", "Illumina HiSeq 4000"]])
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = self.read_rows(out)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["RunAccession"], "SRR1")
        self.assertEqual(rows[0]["Organism"], "Oryza sativa")
        self.assertEqual(rows[0]["SeqType"], "WGS")

    def test_strategy_layers(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "LibraryStrategy"],
                  [["r1", "RNA-Seq"],
                   ["r2", "WGS"],
                   ["r3", "Totally-New-Strategy"]])  # unknown -> Other
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = {r["RunAccession"]: r["SeqType"] for r in self.read_rows(out)}
        self.assertEqual(rows["r1"], "RNA-Seq")
        self.assertEqual(rows["r2"], "WGS")
        self.assertEqual(rows["r3"], "Other")

    def test_source_fallback_genomic_wgs(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "LibraryStrategy", "LibrarySource"],
                  [["r1", "UnknownStrategy", "GENOMIC"]])
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)[0]["SeqType"], "WGS")

    def test_long_read_instrument_layer(self):
        # strategy unknown, no source, PacBio instrument -> Long-Read
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "InstrumentModel"],
                  [["r1", "PacBio Sequel II"]])
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)[0]["SeqType"], "Long-Read")

    def test_title_keyword_layer(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "ExperimentTitle"],
                  [["r1", "Leaf transcriptome rna-seq of rice"]])
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)[0]["SeqType"], "RNA-Seq")

    def test_override_forces_seqtype(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "LibraryStrategy"], [["r1", "WGS"]])
        proc, out = self.run_tool(inp, extra=["--seq-type-override", "RNA-Seq"])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)[0]["SeqType"], "RNA-Seq")

    def test_duplicate_mapping_merged(self):
        # "Run" and "run_accession" both map to RunAccession
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "run_accession"], [["SRR1", "SRR1"]])
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(self.read_rows(out)[0]["RunAccession"], "SRR1")

    def test_ragged_extra_fields_dropped(self):
        # regression: extra fields produced a None key and crashed sorting
        inp = self.write("meta.tsv", "")
        with inp.open("w") as fh:
            fh.write("Run\tLibraryStrategy\n")
            fh.write("SRR1\tWGS\tSURPLUS\tFIELDS\n")  # 2 extra fields
            fh.write("SRR2\tRNA-Seq\n")
        proc, out = self.run_tool(inp)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("extra fields dropped", proc.stderr)
        self.assertEqual(len(self.read_rows(out)), 2)

    def test_stdout_output(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "LibraryStrategy"], [["r1", "WGS"]])
        proc = self.run_script("format_sra_summary.py", "--input", inp,
                               "--output", "-")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("RunAccession", proc.stdout)

    def test_empty_input_exit_1(self):
        inp = self.write("meta.tsv", "")
        write_tsv(inp, ["Run", "LibraryStrategy"], [])
        proc, _ = self.run_tool(inp)
        self.assertEqual(proc.returncode, 1)


if __name__ == "__main__":
    unittest.main()
