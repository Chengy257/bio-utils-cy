"""End-to-end tests for bin/extract_species_rna_rfam.sh.

Uses the real seqkit binary (available on this machine) with tiny gzipped
Rfam-style fixtures; a failing seqkit stub is injected via BUC_SEQKIT_BIN
for the atomicity test.
"""

import gzip
import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import ScriptTestCase  # noqa: E402

SP = "Testus genericus"
OTHER = "Alienus otherus"


class RfamFixtureMixin:

    def make_rfam(self, families):
        """families: {rf_id: (class, [(header, seq), ...])} -> (indir, family.txt)."""
        indir = self.tmp / "rfam"
        indir.mkdir()
        fam_lines = []
        for rf, (cls, seqs) in families.items():
            fam_lines.append(f"{rf}\t{cls}\tdescription of {cls}")
            with gzip.open(indir / f"{rf}.fa.gz", "wt") as fh:
                for i, (hdr, seq) in enumerate(seqs, 1):
                    fh.write(f">{rf}.{i} {hdr}\n{seq}\n")
        family = self.write("family.txt", "\n".join(fam_lines) + "\n")
        return indir, family


class TestExtractSpeciesRnaRfam(ScriptTestCase, RfamFixtureMixin):

    def _run(self, *args):
        return self.run_script("extract_species_rna_rfam.sh", *args)

    def test_basic_extraction_and_empty_class(self):
        indir, family = self.make_rfam({
            "RF00001": ("rRNA", [(f"{SP} embryo", "ACGUACGU"), (OTHER, "UUUU")]),
            "RF00002": ("tRNA", [(OTHER, "GGGG")]),           # species absent
            "RF00003": ("snRNA", [(f"{SP} leaf", "CCCC")]),
        })
        out = self.tmp / "out"
        proc = self._run("-i", str(indir), "-a", str(family), "-o", str(out), SP)
        self.assertEqual(proc.returncode, 0, proc.stderr)

        # all logs on stderr, nothing on stdout
        self.assertEqual(proc.stdout, "")

        # rRNA: only the species' sequence (1 of 2)
        rfa = (out / "Tgen_Rfam_rRNA.fa").read_text()
        self.assertEqual(rfa.count(">"), 1)
        self.assertIn("RF00001.1", rfa)
        self.assertIn(SP, rfa)

        # tRNA: species absent -> seqkit ran, empty file, no crash
        self.assertTrue((out / "Tgen_Rfam_tRNA.fa").exists())
        self.assertEqual((out / "Tgen_Rfam_tRNA.fa").stat().st_size, 0)

        # snRNA extracted; count line reports "1 sequences" (not "0\n0")
        self.assertIn("(1 sequences)", proc.stderr)

        # miRNA class has zero families: grep zero-match must not abort,
        # the .RFids list is empty and no output fasta is produced
        self.assertTrue((out / "miRNA.RFids").exists())
        self.assertEqual((out / "miRNA.RFids").read_text(), "")
        self.assertFalse((out / "Tgen_Rfam_miRNA.fa").exists())
        self.assertIn("No fasta files found for miRNA", proc.stderr)

        # snRNA.RFids contains exactly the snRNA family
        self.assertEqual((out / "snRNA.RFids").read_text(), "RF00003\n")

    def test_species_matched_as_literal_not_regex(self):
        # "Genus. spec" as a regex would also match "GenusX spec other"
        indir, family = self.make_rfam({
            "RF00005": ("rRNA", [("Genus. spec real", "AAAA")]),
            "RF00006": ("rRNA", [("GenusX spec other", "TTTT")]),
        })
        out = self.tmp / "out"
        proc = self._run("-i", str(indir), "-a", str(family),
                         "-o", str(out), "Genus. spec")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rfa = (out / "Gspe_Rfam_rRNA.fa").read_text()
        self.assertEqual(rfa.count(">"), 1)
        self.assertIn("Genus. spec real", rfa)
        self.assertNotIn("GenusX spec other", rfa)

    def test_duplicate_species_processed_once(self):
        indir, family = self.make_rfam({
            "RF00001": ("rRNA", [(f"{SP} x", "ACGU")]),
        })
        out = self.tmp / "out"
        proc = self._run("-i", str(indir), "-a", str(family),
                         "-o", str(out), SP, SP)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stderr.count(f"Processing species: {SP}"), 1)

    def test_short_name_collision_aborts(self):
        indir, family = self.make_rfam({
            "RF00001": ("rRNA", [(SP, "ACGU")]),
        })
        proc = self._run("-i", str(indir), "-a", str(family),
                         "-o", str(self.tmp / "out"),
                         "Testus genericus", "Testus genericus maximus")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("share the short name", proc.stderr)

    def test_invalid_threads_rejected(self):
        indir, family = self.make_rfam({})
        proc = self._run("-i", str(indir), "-a", str(family), "-t", "abc", SP)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("positive integer", proc.stderr)

    def test_seqkit_failure_leaves_no_truncated_output(self):
        indir, family = self.make_rfam({
            "RF00001": ("rRNA", [(SP, "ACGUACGU")]),
        })
        stub = self.tmp / "failing_seqkit.sh"
        stub.write_text("#!/bin/sh\nexit 1\n")
        stub.chmod(0o755)
        out = self.tmp / "out"
        with mock.patch.dict(os.environ, {"BUC_SEQKIT_BIN": str(stub)}):
            proc = self._run("-i", str(indir), "-a", str(family),
                             "-o", str(out), SP)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("seqkit failed", proc.stderr)
        # no output fasta and no temp staging files left behind
        self.assertFalse((out / "Tgen_Rfam_rRNA.fa").exists())
        self.assertEqual(list(out.glob("*.tmp*")), [])

    def test_existing_output_is_skipped(self):
        indir, family = self.make_rfam({
            "RF00001": ("rRNA", [(SP, "ACGU")]),
        })
        out = self.tmp / "out"
        out.mkdir(parents=True)
        pre = out / "Tgen_Rfam_rRNA.fa"
        pre.write_text("PREEXISTING CONTENT\n")
        proc = self._run("-i", str(indir), "-a", str(family),
                         "-o", str(out), SP)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(pre.read_text(), "PREEXISTING CONTENT\n")
        self.assertIn("SKIP", proc.stderr)

    def test_missing_input_dir_fails_cleanly(self):
        family = self.write("family.txt", "RF00001\trRNA\tx\n")
        proc = self._run("-i", str(self.tmp / "nope"), "-a", str(family), SP)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Input directory not found", proc.stderr)


if __name__ == "__main__":
    unittest.main()
