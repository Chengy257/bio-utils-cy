"""Functional tests for bin/extract_utr.py.

Synthetic 200 bp chromosome: 1-50 A, 51-100 C, 101-150 G, 151-200 T.
Expected UTR sequences are hand-derived; negative-strand cases verify
the reverse-complement fix and multi-exon cases verify transcript-order
joining.
"""

import csv
import unittest

from common import ScriptTestCase


def make_genome():
    return ">chr1 test chromosome\n" + "A" * 50 + "C" * 50 + "G" * 50 + "T" * 50 + "\n"


def gtf_line(chrom, feature, start, end, strand, tx_id, extra_attr=""):
    return (f"{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
            f"gene_id \"g_{tx_id}\"; transcript_id \"{tx_id}\";{extra_attr}\n")


class ExtractUtrTest(ScriptTestCase):

    def _run(self, gtf_text, extra_genome=""):
        gtf = self.write("ann.gtf", gtf_text)
        genome = self.write("genome.fa", make_genome() + extra_genome)
        fasta = self.tmp / "utr.fa"
        csvout = self.tmp / "utr.csv"
        proc = self.run_script("extract_utr.py", "--gtf", gtf, "--genome", genome,
                               "-f", fasta, "-c", csvout)
        return proc, fasta, csvout

    @staticmethod
    def _read_fasta(path):
        records = {}
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:]
                records[name] = ""
            else:
                records[name] += line.strip()
        return records

    def test_plus_strand_utr_sequences(self):
        gtf = (gtf_line("chr1", "exon", 10, 60, "+", "t_plus")
               + gtf_line("chr1", "CDS", 30, 60, "+", "t_plus")
               + gtf_line("chr1", "exon", 110, 190, "+", "t_plus")
               + gtf_line("chr1", "CDS", 110, 140, "+", "t_plus"))
        proc, fasta, csvout = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        records = self._read_fasta(fasta)
        self.assertEqual(records["t_plus_5UTR"], "A" * 20)
        self.assertEqual(records["t_plus_3UTR"], "G" * 10 + "T" * 40)

    def test_minus_strand_utr_is_reverse_complemented(self):
        gtf = (gtf_line("chr1", "exon", 10, 80, "-", "t_minus")
               + gtf_line("chr1", "CDS", 30, 60, "-", "t_minus"))
        proc, fasta, _ = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        records = self._read_fasta(fasta)
        # 5'UTR on minus strand = genomic 61-80 (C*20), reverse-complemented.
        self.assertEqual(records["t_minus_5UTR"], "G" * 20)
        # 3'UTR on minus strand = genomic 10-29 (A*20), reverse-complemented.
        self.assertEqual(records["t_minus_3UTR"], "T" * 20)

    def test_multi_exon_utr_joined_in_transcript_order(self):
        # Minus-strand transcript: exons 10-40, 70-100, 130-160; CDS 85-140.
        # 5'UTR = genomic 141-160 (G*10+T*10) -> rc = A*10+C*10.
        # 3'UTR = genomic 70-84 (C*15) + 10-40 (A*31), transcript order is
        # descending genomic on minus strand -> rc(C*15) + rc(A*31).
        gtf = (gtf_line("chr1", "exon", 10, 40, "-", "t_multi")
               + gtf_line("chr1", "exon", 70, 100, "-", "t_multi")
               + gtf_line("chr1", "exon", 130, 160, "-", "t_multi")
               + gtf_line("chr1", "CDS", 85, 100, "-", "t_multi")
               + gtf_line("chr1", "CDS", 130, 140, "-", "t_multi"))
        proc, fasta, _ = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        records = self._read_fasta(fasta)
        self.assertEqual(records["t_multi_5UTR"], "A" * 10 + "C" * 10)
        self.assertEqual(records["t_multi_3UTR"], "G" * 15 + "T" * 31)
        self.assertEqual(len(records), 2)  # one record per type, unique IDs

    def test_csv_lengths(self):
        gtf = (gtf_line("chr1", "exon", 10, 60, "+", "t_plus")
               + gtf_line("chr1", "CDS", 30, 60, "+", "t_plus")
               + gtf_line("chr1", "exon", 110, 190, "+", "t_plus")
               + gtf_line("chr1", "CDS", 110, 140, "+", "t_plus")
               + gtf_line("chr1", "exon", 10, 80, "-", "t_minus")
               + gtf_line("chr1", "CDS", 30, 60, "-", "t_minus")
               + gtf_line("chr1", "exon", 10, 40, "-", "t_multi")
               + gtf_line("chr1", "exon", 70, 100, "-", "t_multi")
               + gtf_line("chr1", "exon", 130, 160, "-", "t_multi")
               + gtf_line("chr1", "CDS", 85, 100, "-", "t_multi")
               + gtf_line("chr1", "CDS", 130, 140, "-", "t_multi"))
        proc, _, csvout = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(csvout) as fh:
            rows = {r["transcript_id"]: r for r in csv.DictReader(fh)}
        self.assertEqual(rows["t_plus"]["5UTR_length"], "20")
        self.assertEqual(rows["t_plus"]["3UTR_length"], "50")
        self.assertEqual(rows["t_minus"]["5UTR_length"], "20")
        self.assertEqual(rows["t_minus"]["3UTR_length"], "20")
        self.assertEqual(rows["t_multi"]["5UTR_length"], "20")
        self.assertEqual(rows["t_multi"]["3UTR_length"], "46")

    def test_noncoding_transcript_skipped(self):
        gtf = (gtf_line("chr1", "exon", 10, 60, "+", "t_plus")
               + gtf_line("chr1", "CDS", 30, 60, "+", "t_plus")
               + gtf_line("chr1", "exon", 10, 90, "+", "t_nc"))  # no CDS
        proc, fasta, csvout = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        records = self._read_fasta(fasta)
        self.assertNotIn("t_nc_5UTR", records)
        with open(csvout) as fh:
            tx_ids = {r["transcript_id"] for r in csv.DictReader(fh)}
        self.assertNotIn("t_nc", tx_ids)
        self.assertIn("t_plus", tx_ids)

    def test_malformed_line_skipped_without_traceback(self):
        bad = "chr1\ttest\texon\tXX\t60\t.\t+\t.\ttranscript_id \"t_bad\";\n"
        gtf = (bad
               + gtf_line("chr1", "exon", 10, 80, "-", "t_minus")
               + gtf_line("chr1", "CDS", 30, 60, "-", "t_minus"))
        proc, fasta, _ = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("malformed", proc.stderr)
        records = self._read_fasta(fasta)
        self.assertEqual(records["t_minus_5UTR"], "G" * 20)

    def test_unquoted_and_versioned_transcript_id(self):
        # Bare (unquoted) transcript_id must parse; transcript_id_version
        # must not shadow transcript_id.
        gtf = (gtf_line("chr1", "exon", 10, 80, "-", "t_bare",
                        extra_attr=" transcript_id_version \"t_bare_v\";")
               + gtf_line("chr1", "CDS", 30, 60, "-", "t_bare",
                          extra_attr=" transcript_id_version \"t_bare_v\";"))
        proc, fasta, _ = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        records = self._read_fasta(fasta)
        self.assertIn("t_bare_5UTR", records)
        self.assertNotIn("t_bare_v_5UTR", records)
        self.assertEqual(records["t_bare_5UTR"], "G" * 20)

    def test_missing_chromosome_warns_and_omits_utr(self):
        gtf = (gtf_line("chr2", "exon", 10, 80, "-", "t_ghost")
               + gtf_line("chr2", "CDS", 30, 60, "-", "t_ghost")
               + gtf_line("chr1", "exon", 10, 80, "-", "t_minus")
               + gtf_line("chr1", "CDS", 30, 60, "-", "t_minus"))
        proc, fasta, csvout = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("not in genome", proc.stderr)
        records = self._read_fasta(fasta)
        self.assertNotIn("t_ghost_5UTR", records)
        self.assertIn("t_minus_5UTR", records)
        with open(csvout) as fh:
            rows = {r["transcript_id"]: r for r in csv.DictReader(fh)}
        self.assertEqual(rows["t_ghost"]["5UTR_length"], "0")


if __name__ == "__main__":
    unittest.main()
