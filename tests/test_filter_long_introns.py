"""Functional tests for bin/filter_long_introns.py.

Pure-Python GTF line filtering replaces the former BCBio.GFF
dependency; the flat-GTF regression (nothing was ever filtered) is a
dedicated test case.
"""

import unittest

from common import ScriptTestCase


def gtf_line(chrom, feature, start, end, strand, gene_id, tx_id):
    return (f"{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
            f"gene_id \"{gene_id}\"; transcript_id \"{tx_id}\";\n")


def gene_line(chrom, start, end, strand, gene_id):
    return (f"{chrom}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
            f"gene_id \"{gene_id}\";\n")


class FilterLongIntronsTest(ScriptTestCase):

    def _run(self, gtf_text, max_intron=20000):
        gtf = self.write("in.gtf", gtf_text)
        out = self.tmp / "filtered.gtf"
        proc = self.run_script("filter_long_introns.py", "-i", gtf, "-o", out,
                               "-m", max_intron)
        return proc, out

    def _standard_gtf(self):
        # g_ok: introns 1000 and 1000 (all <= threshold)
        g_ok = (
            gtf_line("chr1", "transcript", 0, 12000, "+", "g_ok", "t_ok")
            + gtf_line("chr1", "exon", 0, 1000, "+", "g_ok", "t_ok")
            + gtf_line("chr1", "exon", 2000, 3000, "+", "g_ok", "t_ok")
            + gtf_line("chr1", "exon", 11000, 12000, "+", "g_ok", "t_ok")
            + gtf_line("chr1", "CDS", 2000, 3000, "+", "g_ok", "t_ok")
        )
        # g_long: intron 21000-1000 = 20000? no: exon2 starts 31000 -> intron 30000.
        g_long = (
            gtf_line("chr1", "transcript", 0, 32000, "+", "g_long", "t_long")
            + gtf_line("chr1", "exon", 0, 1000, "+", "g_long", "t_long")
            + gtf_line("chr1", "exon", 31000, 32000, "+", "g_long", "t_long")
        )
        header = "# GTF header\n"
        genes = gene_line("chr1", 0, 12000, "+", "g_ok") + gene_line("chr1", 0, 32000, "+", "g_long")
        return header + genes + g_ok + g_long

    def test_long_intron_transcript_removed(self):
        proc, out = self._run(self._standard_gtf(), max_intron=20000)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        text = out.read_text()
        self.assertNotIn("t_long", text)
        self.assertNotIn("g_long", text)      # gene line dropped too
        self.assertIn("t_ok", text)
        self.assertIn("g_ok", text)           # gene line kept
        self.assertIn("# GTF header", text)   # comments kept
        self.assertIn("Total transcripts: 2, Filtered: 1, Retained: 1", proc.stderr)

    def test_retained_lines_are_verbatim(self):
        original = self._standard_gtf()
        proc, out = self._run(original, max_intron=20000)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        out_lines = set(out.read_text().splitlines())
        for line in original.splitlines():
            if "t_long" in line or line.endswith('gene_id "g_long";'):
                self.assertNotIn(line, out_lines)
            else:
                self.assertIn(line, out_lines)  # byte-identical, no reformatting

    def test_boundary_intron_exactly_at_threshold_kept(self):
        # Intron = next.start - prev.end = 21000 - 1000 = 20000 == threshold.
        gtf = (
            gtf_line("chr1", "exon", 0, 1000, "+", "g", "t")
            + gtf_line("chr1", "exon", 21000, 22000, "+", "g", "t")
        )
        proc, out = self._run(gtf, max_intron=20000)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("t", out.read_text())

    def test_intron_one_over_threshold_filtered(self):
        gtf = (
            gtf_line("chr1", "exon", 0, 1000, "+", "g", "t")
            + gtf_line("chr1", "exon", 21001, 22000, "+", "g", "t")
        )
        proc, out = self._run(gtf, max_intron=20000)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "")

    def test_single_exon_transcript_kept(self):
        gtf = gtf_line("chr1", "exon", 0, 50000, "+", "g", "t")
        proc, out = self._run(gtf, max_intron=100)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("t", out.read_text())

    def test_flat_out_of_order_gtf_is_filtered(self):
        # Regression: flat GTFs (exon and gene lines mixed, no hierarchy,
        # exon lines before their gene line) used to bypass the filter.
        gtf = (
            gtf_line("chr1", "exon", 0, 1000, "+", "g_long", "t_long")
            + gtf_line("chr1", "exon", 31000, 32000, "+", "g_long", "t_long")
            + gene_line("chr1", 0, 32000, "+", "g_long")
            + gtf_line("chr1", "exon", 0, 5000, "+", "g_ok", "t_ok")
            + gtf_line("chr1", "exon", 6000, 9000, "+", "g_ok", "t_ok")
        )
        proc, out = self._run(gtf, max_intron=20000)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        text = out.read_text()
        self.assertNotIn("t_long", text)
        self.assertNotIn("g_long", text)
        self.assertIn("t_ok", text)

    def test_no_transcript_id_exits_1(self):
        gtf = gene_line("chr1", 0, 1000, "+", "g1")
        proc, _ = self._run(gtf, max_intron=20000)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("no transcript_id", proc.stderr)

    def test_empty_input_exits_1(self):
        proc, _ = self._run("")
        self.assertEqual(proc.returncode, 1)

    def test_all_transcripts_removed_yields_empty_output(self):
        gtf = (
            gtf_line("chr1", "exon", 0, 1000, "+", "g", "t1")
            + gtf_line("chr1", "exon", 99000, 100000, "+", "g", "t1")
        )
        proc, out = self._run(gtf, max_intron=20000)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out.read_text(), "")


if __name__ == "__main__":
    unittest.main()
