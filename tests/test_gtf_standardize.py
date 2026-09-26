"""Functional tests for bin/gtf_standardize.sh.

Covers both processing modes (with/without gene lines), the
out-of-order backfill regression, natural chromosome sorting, and the
POSIX-awk attribute parsing.
"""

import unittest

from common import ScriptTestCase


def line(chrom, feat, start, end, strand, attrs):
    return f"{chrom}\ttest\t{feat}\t{start}\t{end}\t.\t{strand}\t0\t{attrs}\n"


def parse(path):
    rows = {}
    for ln in path.read_text().splitlines():
        f = ln.split("\t")
        rows.setdefault(f[2], []).append(f)
    return rows


class GtfStandardizeTest(ScriptTestCase):

    def _run(self, gtf_text, *extra):
        gtf = self.write("in.gtf", gtf_text)
        out = self.tmp / "out.gtf"
        proc = self.run_script("gtf_standardize.sh", "-i", gtf, "-o", out, *extra)
        return proc, out

    def test_with_gene_mode(self):
        gtf = (
            line("chr1", "gene", 100, 900, "+", 'gene_id "g1"; gene_name "G One"; gene_biotype "protein_coding";')
            + line("chr1", "transcript", 100, 900, "+", 'gene_id "g1"; transcript_id "t1"; gene_biotype "protein_coding";')
            + line("chr1", "exon", 100, 300, "+", 'gene_id "g1"; transcript_id "t1"; exon_number "1"; gene_biotype "protein_coding";')
            + line("chr1", "CDS", 200, 300, "+", 'gene_id "g1"; transcript_id "t1"; exon_number "1"; protein_id "p1"; gene_biotype "protein_coding";')
            + line("chr1", "Selenocysteine", 250, 252, "+", 'gene_id "g1"; transcript_id "t1";')
        )
        proc, out = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = parse(out)
        self.assertNotIn("Selenocysteine", rows)  # unsupported type dropped
        gene = rows["gene"][0]
        self.assertEqual(gene[8], 'gene_id "g1"; gene_name "G One"; gene_biotype "protein_coding"')
        tx = rows["transcript"][0]
        # gene_name backfilled from the gene line
        self.assertIn('gene_name "G One"', tx[8])
        cds = rows["CDS"][0]
        self.assertEqual(cds[7], "0")  # phase preserved
        self.assertIn('protein_id "p1"', cds[8])
        exon = rows["exon"][0]
        self.assertNotIn("protein_id", exon[8])

    def test_without_gene_mode_synthesizes_and_backfills_out_of_order(self):
        # Exon/CDS lines come BEFORE their transcript lines: the former
        # implementation silently lost gene_id on such lines.
        gtf = (
            line("chr1", "exon", 500, 900, "+", 'transcript_id "t2"; exon_number "2";')
            + line("chr1", "CDS", 200, 300, "+", 'transcript_id "t1"; exon_number "1"; protein_id "p1";')
            + line("chr1", "transcript", 100, 900, "+", 'gene_id "gX"; transcript_id "t1"; gene_biotype "protein_coding"; gene_name "GX";')
            + line("chr1", "transcript", 500, 900, "+", 'gene_id "gX"; transcript_id "t2"; gene_biotype "protein_coding";')
            + line("chr1", "exon", 100, 300, "+", 'transcript_id "t1"; exon_number "1";')
        )
        proc, out = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = parse(out)
        self.assertIn("gene", rows)
        gene = rows["gene"][0]
        # gene extent = min/max over its transcript lines
        self.assertEqual([gene[3], gene[4]], ["100", "900"])
        self.assertIn('gene_biotype "protein_coding"', gene[8])
        self.assertIn('gene_name "GX"', gene[8])
        for feat in ("exon", "CDS"):
            for r in rows[feat]:
                self.assertIn('gene_id "gX"', r[8])
        cds = rows["CDS"][0]
        self.assertIn('gene_biotype "protein_coding"', cds[8])  # backfilled
        self.assertIn('gene_name "GX"', cds[8])

    def test_natural_chromosome_sorting(self):
        gtf = (
            line("chr10", "exon", 10, 20, "+", 'gene_id "ga"; transcript_id "ta"; gene_biotype "b";')
            + line("chr2", "exon", 10, 20, "+", 'gene_id "gb"; transcript_id "tb"; gene_biotype "b";')
            + line("chr1", "exon", 10, 20, "+", 'gene_id "gc"; transcript_id "tc"; gene_biotype "b";')
        )
        proc, out = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        chroms = [ln.split("\t")[0] for ln in out.read_text().splitlines()]
        self.assertEqual(chroms, ["chr1", "chr2", "chr10"])  # not chr1, chr10, chr2

    def test_custom_attrs(self):
        gtf = line("chr1", "exon", 10, 20, "+", 'gene_id "g1"; transcript_id "t1"; gene_biotype "b"; gene_name "N";')
        proc, out = self._run(gtf, "-a", "gene_id,transcript_id")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        attrs = parse(out)["exon"][0][8]
        self.assertEqual(attrs, 'gene_id "g1"; transcript_id "t1"')

    def test_unknown_attr_rejected(self):
        gtf = line("chr1", "exon", 10, 20, "+", 'gene_id "g1";')
        proc, _ = self._run(gtf, "-a", "gene_id,chromosome")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Unknown attribute", proc.stderr)

    def test_overwrite_protection(self):
        gtf = line("chr1", "exon", 10, 20, "+", 'gene_id "g1"; transcript_id "t1"; gene_biotype "b";')
        proc, out = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        proc, _ = self._run(gtf)  # same -o, no -f
        self.assertEqual(proc.returncode, 1)
        self.assertIn("already exists", proc.stderr)
        proc, _ = self._run(gtf, "-f")
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_logs_go_to_stderr(self):
        gtf = line("chr1", "exon", 10, 20, "+", 'gene_id "g1"; transcript_id "t1"; gene_biotype "b";')
        proc, _ = self._run(gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, "")

    def test_missing_input(self):
        proc = self.run_script("gtf_standardize.sh", "-i", "/nonexistent.gtf")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("does not exist", proc.stderr)

    def test_default_output_suffix(self):
        gtf = self.write("ann.gtf", line("chr1", "exon", 10, 20, "+", 'gene_id "g1"; transcript_id "t1"; gene_biotype "b";'))
        proc = self.run_script("gtf_standardize.sh", "-i", gtf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertTrue((self.tmp / "ann.gtf.standardized.gtf").exists())


if __name__ == "__main__":
    unittest.main()
