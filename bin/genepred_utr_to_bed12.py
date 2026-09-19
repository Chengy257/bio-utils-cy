#!/usr/bin/env python3
#########################################################################
# File Name: genepred_utr_to_bed12.py
# Author: ChengYu
# Description: Extract UTR regions from genePred annotations and
#              output as BED12 format.
# Created Time: 2026
#########################################################################
"""Extract UTR regions from genePred annotations as BED12.

Converts genePred format annotations to BED12 files containing only
the UTR portions (5'UTR, 3'UTR, or both). Supports multi-threaded
processing.
"""

import argparse
import logging
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import List, Optional, Tuple

__version__ = "1.0.0"


def parse_genepred_line(line: str) -> Tuple:
    """Parse a single genePred line.

    Args:
        line: Tab-separated genePred line.

    Returns:
        Tuple of (gene_id, chrom, strand, tx_start, tx_end, cds_start,
                  cds_end, exon_starts, exon_ends).
    """
    fields = line.strip().split("\t")
    return (
        fields[0],                          # gene_id
        fields[1],                          # chrom
        fields[2],                          # strand
        int(fields[3]),                     # tx_start
        int(fields[4]),                     # tx_end
        int(fields[5]),                     # cds_start
        int(fields[6]),                     # cds_end
        list(map(int, fields[8].rstrip(",").split(","))),  # exon_starts
        list(map(int, fields[9].rstrip(",").split(","))),  # exon_ends
    )


def extract_utr_bed12(
    gene_id: str,
    chrom: str,
    strand: str,
    tx_start: int,
    tx_end: int,
    cds_start: int,
    cds_end: int,
    exon_starts: List[int],
    exon_ends: List[int],
    utr_type: str,
) -> Optional[List]:
    """Extract UTR region and format as BED12.

    Args:
        gene_id through exon_ends: genePred fields.
        utr_type: '5UTR' or '3UTR'.

    Returns:
        BED12 row as list, or None if no UTR for this type.
    """
    block_sizes = []
    block_starts = []
    chrom_start = None
    chrom_end = None

    if utr_type == "5UTR":
        for es, ee in zip(exon_starts, exon_ends):
            if ee <= cds_start:
                if chrom_start is None:
                    chrom_start = es
                chrom_end = ee
                block_sizes.append(ee - es)
                block_starts.append(es - chrom_start)
            elif es < cds_start:
                if chrom_start is None:
                    chrom_start = es
                chrom_end = cds_start
                block_sizes.append(cds_start - es)
                block_starts.append(es - chrom_start)
                break
            else:
                break
    elif utr_type == "3UTR":
        for es, ee in zip(exon_starts, exon_ends):
            if es >= cds_end:
                if chrom_start is None:
                    chrom_start = es
                chrom_end = ee
                block_sizes.append(ee - es)
                block_starts.append(es - chrom_start)
            elif ee > cds_end:
                if chrom_start is None:
                    chrom_start = cds_end
                chrom_end = ee
                block_sizes.append(ee - cds_end)
                block_starts.append(cds_end - chrom_start)

    if chrom_start is None or not block_sizes:
        return None

    return [
        chrom, chrom_start, chrom_end, gene_id, 0, strand,
        chrom_start, chrom_start, 0,
        len(block_sizes),
        ",".join(map(str, block_sizes)) + ",",
        ",".join(map(str, block_starts)) + ",",
    ]


def process_line(args: tuple) -> Optional[List]:
    """Process a single genePred line for a given UTR type."""
    line, utr_type = args
    try:
        gene_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends = parse_genepred_line(line)
        return extract_utr_bed12(
            gene_id, chrom, strand, tx_start, tx_end,
            cds_start, cds_end, exon_starts, exon_ends, utr_type,
        )
    except (IndexError, ValueError) as e:
        logging.warning("Failed to parse line: %s", e)
        return None


def convert(input_file: str, output_file: str, utr_type: str, threads: int = 1) -> None:
    """Convert genePred to BED12 with UTR regions.

    Args:
        input_file: Input genePred file.
        output_file: Output BED12 file.
        utr_type: '5UTR', '3UTR', or 'both'.
        threads: Number of parallel workers.
    """
    with open(input_file, "r") as fh:
        lines = [line for line in fh if line.strip()]

    if not lines:
        logging.error("Empty input file: %s", input_file)
        sys.exit(1)

    types = ["5UTR", "3UTR"] if utr_type == "both" else [utr_type]

    for ut in types:
        if utr_type == "both":
            out_path = output_file.replace(".bed12", f"_{ut}.bed12").replace(".bed", f"_{ut}.bed")
            if out_path == output_file:
                out_path = f"{output_file}.{ut}"
        else:
            out_path = output_file

        if threads > 1:
            with Pool(threads) as pool:
                results = pool.map(process_line, [(line, ut) for line in lines])
        else:
            results = [process_line((line, ut)) for line in lines]

        count = 0
        with open(out_path, "w") as fh:
            for row in results:
                if row:
                    fh.write("\t".join(map(str, row)) + "\n")
                    count += 1

        logging.info("%s: %d regions -> %s", ut, count, out_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract UTR regions from genePred annotations as BED12 format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python genepred_utr_to_bed12.py -i annotation.gp -o utr.bed12 --utr 5UTR
  python genepred_utr_to_bed12.py -i annotation.gp -o utr.bed12 --utr both -t 4
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input genePred file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output BED12 file.")
    parser.add_argument(
        "--utr", type=str, choices=["5UTR", "3UTR", "both"], default="5UTR",
        help="UTR type to extract (default: 5UTR).",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1).")
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

    if not Path(args.input).is_file():
        logging.error("Input file not found: %s", args.input)
        sys.exit(1)

    convert(args.input, args.output, args.utr, args.threads)


if __name__ == "__main__":
    main()
