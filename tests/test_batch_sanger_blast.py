"""Functional tests for bin/batch_sanger_blast.sh.

Runs the real blastn/makeblastdb (skipped when BLAST+ is not on PATH)
against small synthetic Sanger exports.
"""

import os
import shutil
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ScriptTestCase  # noqa: E402

BLASTN = shutil.which("blastn")
MAKEBLASTDB = shutil.which("makeblastdb")

_REF_SEG = ("ACGTTGCAAGCTTACGGATCGTAGCTAGCTTACGGATCCATGGACTGACCTTACGGATCCATG"
            "GACTGACCTTAGGCATGCCATTACGGATCCATGGACTGACCTTACGGATCCATGGACTGAC")
REF_500 = (_REF_SEG * 5)[:500]
assert len(REF_500) == 500


@unittest.skipUnless(BLASTN and MAKEBLASTDB, "BLAST+ not available")
class TestBatchSangerBlast(ScriptTestCase):
    def _write_ref(self):
        return self.write("ref.fa", ">ref_chr1\n%s\n" % REF_500)

    def _workdir(self, files):
        """files: {relpath: content} inside a 'sanger_zips' dir."""
        work = self.tmp / "sanger_zips"
        work.mkdir(parents=True, exist_ok=True)
        for rel, content in files.items():
            p = work / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(content)
        return work

    def test_full_run_outputs_and_hit(self):
        ref = self._write_ref()
        work = self._workdir({"A1.seq": REF_500[100:400] + "\n"})
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work,
                               "-o", "sanger_BLAST")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        tsv = work / "sanger_BLAST_result.tsv"
        txt = work / "sanger_BLAST_result.txt"
        self.assertTrue(tsv.exists() and txt.exists(), proc.stderr)
        cols = tsv.read_text().split("\t")
        self.assertEqual(cols[0], "A1")
        self.assertEqual(cols[1], "ref_chr1")
        self.assertEqual(cols[2], "100.000")

    def test_duplicate_basenames_get_unique_queries(self):
        ref = self._write_ref()
        work = self._workdir({
            "subA/A1.seq": REF_500[100:400] + "\n",
            "subB/A1.seq": REF_500[0:250] + "\n",
        })
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        fasta = (work / "allseq.fa").read_text()
        self.assertIn(">subA_A1", fasta)
        self.assertIn(">subB_A1", fasta)
        self.assertNotIn(">A1\n", fasta)
        tsv = (work / "sanger_BLAST_result.tsv").read_text().splitlines()
        queries = {ln.split("\t")[0] for ln in tsv if ln.strip()}
        self.assertEqual(queries, {"subA_A1", "subB_A1"})

    def test_non_sequence_file_skipped(self):
        ref = self._write_ref()
        work = self._workdir({
            "good.seq": REF_500[100:400] + "\n",
            "notes.myseq": "lab notes, definitely not a sequence\n",
        })
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("not recognisable as nucleotide sequence", proc.stderr)
        fasta = (work / "allseq.fa").read_text()
        self.assertIn(">good", fasta)
        self.assertNotIn("notes", fasta)

    def test_embedded_fasta_header_stripped(self):
        ref = self._write_ref()
        work = self._workdir({
            "fasta.seq": ">embedded_header\n%s\n" % REF_500[100:400],
        })
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        fasta = (work / "allseq.fa").read_text()
        self.assertEqual(fasta.count(">"), 1)  # only our own query header
        self.assertIn(">fasta", fasta)
        tsv = (work / "sanger_BLAST_result.tsv").read_text()
        self.assertIn("fasta\tref_chr1", tsv)

    def test_corrupt_zip_warns_but_does_not_abort(self):
        ref = self._write_ref()
        work = self._workdir({
            "broken.zip": "",
            "A1.seq": REF_500[100:400] + "\n",
        })
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Failed to unzip", proc.stderr)
        self.assertTrue((work / "sanger_BLAST_result.tsv").exists())

    def test_real_zip_is_extracted(self):
        if not shutil.which("zip"):
            self.skipTest("zip not available")
        import zipfile
        ref = self._write_ref()
        work = self._workdir({})
        with zipfile.ZipFile(work / "pkg.zip", "w") as zf:
            zf.writestr("inside.seq", REF_500[350:500] + "\n")
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        fasta = (work / "allseq.fa").read_text()
        self.assertIn(">inside", fasta)
        self.assertIn("Extracted 1 archive", proc.stderr)

    def test_relative_db_path_resolved_before_cd(self):
        # documented usage "-d ref.fa -w sanger_zips/" from the parent dir
        ref = self._write_ref()
        work = self._workdir({"A1.seq": REF_500[100:400] + "\n"})
        old_cwd = os.getcwd()
        os.chdir(self.tmp)
        try:
            proc = self.run_script("batch_sanger_blast.sh", "-d", "ref.fa",
                                   "-w", "sanger_zips")
        finally:
            os.chdir(old_cwd)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((work / "sanger_BLAST_result.tsv").exists())

    def test_no_valid_seq_files_fails(self):
        ref = self._write_ref()
        work = self._workdir({})
        proc = self.run_script("batch_sanger_blast.sh", "-d", ref, "-w", work)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No valid *seq files", proc.stderr)

    def test_missing_required_args(self):
        ref = self._write_ref()
        work = self._workdir({"A1.seq": REF_500[100:400] + "\n"})

        proc = self.run_script("batch_sanger_blast.sh", "-h")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("Usage:", proc.stdout)

        proc = self.run_script("batch_sanger_blast.sh", "-w", work)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("-d", proc.stderr)

        proc = self.run_script("batch_sanger_blast.sh", "-d", ref)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("-w", proc.stderr)

        proc = self.run_script("batch_sanger_blast.sh", "-d", self.tmp / "nope.fa",
                               "-w", work)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Database file not found", proc.stderr)

    def test_invalid_threads_rejected(self):
        ref = self._write_ref()
        work = self._workdir({"A1.seq": REF_500[100:400] + "\n"})
        for bad in ("0", "-3", "abc"):
            proc = self.run_script("batch_sanger_blast.sh", "-d", ref,
                                   "-w", work, "-t", bad)
            self.assertEqual(proc.returncode, 1, bad)
            self.assertIn("positive integer", proc.stderr)


if __name__ == "__main__":
    unittest.main()
