#!/usr/bin/env python3
#########################################################################
# File Name: maf_extract_regions.py
# Author: ChengYu
# Description: Extract genomic regions from MAF alignment files based
#              on BED12 coordinates and translate to protein.
# Created Time: 2026
#########################################################################
"""Extract genomic regions from MAF alignment files using BED12 coordinates.

For each region in a BED12 file, extracts the aligned sequence from a
MAF file for a reference species, optionally reverse-complementing for
negative strand. Outputs both nucleotide and translated protein FASTA.

Requires: biopython, bx-python
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__version__ = "1.0.0"


def read_bed12(bed_file: str) -> Dict[str, dict]:
    """Parse a BED12 file into region dictionaries.

    Args:
        bed_file: Path to BED12 file.

    Returns:
        Dict mapping region name -> {chrom, strand, blocks}.
    """
    regions = {}
    with open(bed_file) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < 12:
                logging.warning("Skipping incomplete BED12 line: %s", line.strip())
                continue

            chrom = fields[0]
            start = int(fields[1])
            name = fields[3]
            strand = fields[5]
            block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]

            regions[name] = {
                "chrom": chrom,
                "strand": strand,
                "blocks": [
                    (start + bs, start + bs + sz)
                    for bs, sz in zip(block_starts, block_sizes)
                ],
            }
    logging.info("Loaded %d regions from %s", len(regions), bed_file)
    return regions


def extract_sequences(
    maf_file: str,
    regions: Dict[str, dict],
    reference_species: str,
) -> Dict[str, str]:
    """Extract sequences from MAF for each BED12 region.

    Args:
        maf_file: Path to MAF alignment file.
        regions: BED12 region dicts from read_bed12().
        reference_species: Reference species prefix (e.g., 'hg38').

    Returns:
        Dict mapping region name -> nucleotide sequence string.
    """
    from bx.align import maf

    sequences = {}

    with open(maf_file) as fh:
        reader = maf.Reader(fh)

        for alignment in reader:
            for name, info in regions.items():
                if name in sequences:
                    continue

                ref_prefix = f"{reference_species}.{info['chrom']}"

                for comp in alignment.components:
                    if not comp.src.startswith(ref_prefix):
                        continue

                    cds_seq = ""
                    try:
                        src_parts = comp.src.split(":")
                        if len(src_parts) < 2:
                            continue
                        coords = src_parts[1]
                        aln_start, aln_end = map(int, coords.split("-"))

                        for block_start, block_end in info["blocks"]:
                            # Check overlap
                            if block_end <= aln_start or block_start >= aln_end:
                                continue

                            ov_start = max(block_start, aln_start)
                            ov_end = min(block_end, aln_end)
                            seq_start = ov_start - aln_start
                            seq_end = ov_end - aln_start
                            cds_seq += comp.text[seq_start:seq_end]

                    except (ValueError, IndexError) as e:
                        logging.debug("Error parsing component %s: %s", comp.src, e)
                        continue

                    if cds_seq:
                        if info["strand"] == "-":
                            cds_seq = str(Seq(cds_seq).reverse_complement())
                        sequences[name] = cds_seq.replace("-", "")
                    break  # Found matching component

    logging.info("Extracted %d/%d regions.", len(sequences), len(regions))
    return sequences


def write_fasta(sequences: Dict[str, str], output_path: str, translate: bool = False) -> None:
    """Write sequences to FASTA file.

    Args:
        sequences: Dict of name -> sequence.
        output_path: Output FASTA path.
        translate: If True, translate nucleotide to protein.
    """
    records = []
    for name, seq in sequences.items():
        if translate:
            seq = str(Seq(seq).translate(to_stop=True))
        records.append(SeqRecord(Seq(seq), id=name, description=""))
    SeqIO.write(records, output_path, "fasta")
    logging.info("Wrote %d sequences -> %s", len(records), output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract genomic regions from MAF using BED12 coordinates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python maf_extract_regions.py --bed cds.bed12 --maf alignment.maf --ref hg38 --cds output.fa --protein protein.fa
""",
    )
    parser.add_argument("--bed", type=str, required=True, help="Input BED12 file with regions to extract.")
    parser.add_argument("--maf", type=str, required=True, help="Input MAF alignment file.")
    parser.add_argument("--ref", type=str, required=True, help="Reference species name (e.g., hg38, mm10).")
    parser.add_argument("--cds", type=str, required=True, help="Output nucleotide FASTA file.")
    parser.add_argument("--protein", type=str, required=True, help="Output protein FASTA file.")
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=getattr(logging, args.log_level),
    )

    for label, path in [("BED12", args.bed), ("MAF", args.maf)]:
        if not Path(path).is_file():
            logging.error("%s file not found: %s", label, path)
            sys.exit(1)

    try:
        from bx.align import maf  # noqa: F401
    except ImportError:
        logging.error("Missing dependency: bx-python. Install with: pip install bx-python")
        sys.exit(1)

    regions = read_bed12(args.bed)
    sequences = extract_sequences(args.maf, regions, args.ref)
    write_fasta(sequences, args.cds, translate=False)
    write_fasta(sequences, args.protein, translate=True)


if __name__ == "__main__":
    main()
