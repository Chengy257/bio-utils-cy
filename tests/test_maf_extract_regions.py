"""Functional tests for bin/maf_extract_regions.py.

Mini MAF fixtures with hand-derived expected sequences. These exercise
the four former correctness bugs: gap-column misalignment, cross-block
truncation, minus-strand components, and src prefix collisions — plus
the replaced bx-python parsing (s-line start/size fields).
"""

import unittest

from common import ScriptTestCase


def maf_block(components):
    """components: list of (src, start, size, strand, text)."""
    lines = ["a score=0\n"]
    for src, start, size, strand, text in components:
        lines.append(f"s {src} {start} {size} {strand} 1000 {text}\n")
    lines.append("\n")
    return "".join(lines)


def bed12_line(chrom, start, end, name, strand, blocks):
    """blocks: absolute (bstart, bend) intervals."""
    sizes = [e - s for s, e in blocks]
    starts = [s - start for s, e in blocks]
    return (f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t0\t"
            f"{len(sizes)}\t" + ",".join(map(str, sizes)) + ",\t"
            + ",".join(map(str, starts)) + ",\n")


def read_fasta(path):
    records = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


class MafExtractRegionsTest(ScriptTestCase):

    def _run(self, bed_text, maf_text):
        bed = self.write("regions.bed12", bed_text)
        maf = self.write("aln.maf", maf_text)
        cds = self.tmp / "cds.fa"
        prot = self.tmp / "prot.fa"
        proc = self.run_script("maf_extract_regions.py", "--bed", bed, "--maf", maf,
                               "--ref", "hg38", "--cds", cds, "--protein", prot)
        return proc, cds, prot

    def test_basic_plus_strand_and_translation(self):
        bed = bed12_line("chr1", 10, 20, "r1", "+", [(10, 20)])
        maf = maf_block([
            ("hg38.chr1", 10, 10, "+", "ACGTACGTAC"),
            ("mm10.chr1", 10, 10, "+", "ACGTTCGTAC"),
        ])
        proc, cds, prot = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["r1"], "ACGTACGTAC")
        # ACG TAC GTA -> T Y V (trailing partial codon dropped).
        self.assertEqual(read_fasta(prot)["r1"], "TYV")

    def test_reference_gap_columns_do_not_misalign(self):
        # The reference carries two gap columns; positions 10-17 map to
        # A C G T A C G T, skipping gaps.
        bed = bed12_line("chr1", 10, 18, "r1", "+", [(10, 18)])
        maf = maf_block([
            ("hg38.chr1", 10, 8, "+", "ACGT-A-CGT"),
            ("mm10.chr1", 10, 10, "+", "ACGTAAACGT"),
        ])
        proc, cds, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["r1"], "ACGTACGT")

    def test_region_spanning_multiple_blocks(self):
        bed = bed12_line("chr1", 10, 35, "r1", "+", [(10, 15), (30, 35)])
        maf = (maf_block([("hg38.chr1", 10, 5, "+", "ACGTA")])
               + maf_block([("hg38.chr1", 30, 5, "+", "TTTAA")]))
        proc, cds, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["r1"], "ACGTATTTAA")

    def test_minus_strand_region_reverse_complemented(self):
        bed = bed12_line("chr1", 10, 20, "r1", "-", [(10, 20)])
        maf = maf_block([("hg38.chr1", 10, 10, "+", "ACGTACGTAC")])
        proc, cds, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["r1"], "GTACGTACGT")

    def test_minus_strand_reference_component(self):
        # Plus-strand genome [10,20) is ACGTACGTAC; a minus-strand
        # component stores its reverse complement as text.
        bed_plus = bed12_line("chr1", 10, 20, "rp", "+", [(10, 20)])
        bed_minus = bed12_line("chr1", 10, 20, "rm", "-", [(10, 20)])
        maf = maf_block([("hg38.chr1", 10, 10, "-", "GTACGTACGT")])
        proc, cds, _ = self._run(bed_plus, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["rp"], "ACGTACGTAC")
        proc, cds, _ = self._run(bed_minus, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["rm"], "GTACGTACGT")

    def test_no_chr1_chr10_prefix_collision(self):
        bed = bed12_line("chr1", 10, 20, "r1", "+", [(10, 20)])
        maf = (maf_block([("hg38.chr10", 0, 10, "+", "GGGGGGGGGG")])
               + maf_block([("hg38.chr1", 10, 10, "+", "ACGTACGTAC")]))
        proc, cds, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["r1"], "ACGTACGTAC")

    def test_uncovered_positions_filled_with_n(self):
        bed = bed12_line("chr1", 5, 20, "r1", "+", [(5, 20)])
        maf = maf_block([("hg38.chr1", 10, 10, "+", "ACGTACGTAC")])
        proc, cds, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(read_fasta(cds)["r1"], "NNNNNACGTACGTAC")

    def test_wrong_ref_species_exits_1(self):
        bed = bed12_line("chr1", 10, 20, "r1", "+", [(10, 20)])
        maf = maf_block([("mm10.chr1", 10, 10, "+", "ACGTACGTAC")])
        proc, _, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("No regions found", proc.stderr)

    def test_duplicate_bed_name_keeps_last(self):
        bed = (bed12_line("chr1", 10, 20, "dup", "+", [(10, 20)])
               + bed12_line("chr1", 50, 60, "dup", "+", [(50, 60)]))
        maf = (maf_block([("hg38.chr1", 10, 10, "+", "ACGTACGTAC")])
               + maf_block([("hg38.chr1", 50, 10, "+", "TTTTTTTTTT")]))
        proc, cds, _ = self._run(bed, maf)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("Duplicate region name", proc.stderr)
        records = read_fasta(cds)
        self.assertEqual(list(records), ["dup"])
        self.assertEqual(records["dup"], "TTTTTTTTTT")


if __name__ == "__main__":
    unittest.main()
