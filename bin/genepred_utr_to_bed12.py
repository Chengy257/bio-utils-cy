#!/usr/bin/env python3
#########################################################################
# File Name: genepred_utr_to_bed12.py
# Author: ChengYu
# Description: Extract UTR regions from genePred annotations and
#              output as BED12 format.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: negative-strand 5'/3' UTRs were swapped (interval selection
#     used plus-strand logic for both strands); UTR sides are now chosen
#     by strand: on '-' the 5'UTR lies right of cdsEnd and the 3'UTR
#     left of cdsStart.
#   - FIX: bin-prefixed 11-column genePred (UCSC knownGene exports) was
#     read with an 8-column offset, silently producing garbage; the
#     leading bin column is now auto-detected and stripped, and exon
#     arrays are validated against exonCount.
#   - FIX: non-coding rows (cdsStart == cdsEnd) are skipped with a
#     warning instead of emitting the whole flanking region as "UTR".
#   - FIX: --utr both no longer produces a double suffix for output
#     paths ending in .bed12.
#   - Malformed lines are counted; the tool exits 1 when no line parses.
#########################################################################
"""Extract UTR regions from genePred annotations as BED12.

Converts genePred format annotations to BED12 files containing only
the UTR portions (5'UTR, 3'UTR, or both). Supports multi-threaded
processing.

UTR side by strand:
  plus strand : 5'UTR = left of cdsStart, 3'UTR = right of cdsEnd
  minus strand: 5'UTR = right of cdsEnd,  3'UTR = left of cdsStart

Non-coding rows (cdsStart == cdsEnd) are skipped. Both the standard
10-column genePred and bin-prefixed 11-column exports (UCSC knownGene)
are accepted.
"""

import argparse
import logging
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import List, Optional, Tuple

__version__ = "1.1.0"


def parse_genepred_line(line: str) -> Tuple:
    """Parse a single genePred line (10-column standard, or bin-prefixed
    11-column UCSC knownGene export).

    Args:
        line: Tab-separated genePred line.

    Returns:
        Tuple of (gene_id, chrom, strand, tx_start, tx_end, cds_start,
                  cds_end, exon_starts, exon_ends).

    Raises:
        ValueError: if the line has too few columns or inconsistent
            exon arrays.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) >= 11 and fields[0].isdigit():
        fields = fields[1:]  # strip leading bin column
    if len(fields) < 10:
        raise ValueError(f"expected >=10 genePred columns, got {len(fields)}")

    exon_starts = [int(x) for x in fields[8].rstrip(",").split(",") if x]
    exon_ends = [int(x) for x in fields[9].rstrip(",").split(",") if x]
    if not exon_starts or len(exon_starts) != len(exon_ends):
        raise ValueError(f"inconsistent exon_starts/exon_ends for {fields[0]}")
    if fields[7].isdigit() and int(fields[7]) != len(exon_starts):
        raise ValueError(f"exonCount != len(exon_starts) for {fields[0]}")

    return (
        fields[0],                          # gene_id
        fields[1],                          # chrom
        fields[2],                          # strand
        int(fields[3]),                     # tx_start
        int(fields[4]),                     # tx_end
        int(fields[5]),                     # cds_start
        int(fields[6]),                     # cds_end
        exon_starts,
        exon_ends,
    )


def _utr_blocks(
    exon_starts: List[int], exon_ends: List[int], cut: int, side: str
) -> List[Tuple[int, int]]:
    """Exon portions strictly on one side of a CDS boundary (0-based
    half-open genomic coordinates, ascending order).

    side='left' : portions with end <= cut, plus the partial exon cut at cut.
    side='right': portions with start >= cut, plus the partial exon cut at cut.
    """
    blocks = []
    if side == "left":
        for es, ee in zip(exon_starts, exon_ends):
            if ee <= cut:
                blocks.append((es, ee))
            elif es < cut:
                blocks.append((es, cut))
                break
            else:
                break
    else:
        for es, ee in zip(exon_starts, exon_ends):
            if es >= cut:
                blocks.append((es, ee))
            elif ee > cut:
                blocks.append((cut, ee))
    return blocks


def extract_utr_bed12(
    gene_id: str,
    chrom: str,
    strand: str,
    cds_start: int,
    cds_end: int,
    exon_starts: List[int],
    exon_ends: List[int],
    utr_type: str,
) -> Optional[List]:
    """Extract UTR region and format as BED12.

    Args:
        gene_id, chrom, strand, cds_start, cds_end, exon_starts,
            exon_ends: genePred fields.
        utr_type: '5UTR' or '3UTR'.

    Returns:
        BED12 row as list, or None if no UTR for this type.
    """
    if utr_type == "5UTR":
        blocks = (
            _utr_blocks(exon_starts, exon_ends, cds_end, "right")
            if strand == "-"
            else _utr_blocks(exon_starts, exon_ends, cds_start, "left")
        )
    else:
        blocks = (
            _utr_blocks(exon_starts, exon_ends, cds_start, "left")
            if strand == "-"
            else _utr_blocks(exon_starts, exon_ends, cds_end, "right")
        )

    if not blocks:
        return None

    chrom_start = blocks[0][0]
    chrom_end = blocks[-1][1]
    block_sizes = [ee - es for es, ee in blocks]
    block_starts = [es - chrom_start for es, ee in blocks]

    return [
        chrom, chrom_start, chrom_end, gene_id, 0, strand,
        chrom_start, chrom_start, 0,
        len(blocks),
        ",".join(map(str, block_sizes)) + ",",
        ",".join(map(str, block_starts)) + ",",
    ]


def process_line(args: tuple) -> Tuple[str, Optional[List]]:
    """Process a single genePred line for a given UTR type.

    Returns:
        (status, row) where status is one of 'ok', 'no_utr',
        'noncoding', or 'parse_error'.
    """
    line, utr_type = args
    try:
        (gene_id, chrom, strand, tx_start, tx_end,
         cds_start, cds_end, exon_starts, exon_ends) = parse_genepred_line(line)
    except (IndexError, ValueError) as e:
        logging.warning("Failed to parse line: %s", e)
        return ("parse_error", None)

    if cds_start == cds_end:
        logging.warning("Skipping non-coding transcript %s (cdsStart == cdsEnd).", gene_id)
        return ("noncoding", None)
    if not (tx_start <= cds_start <= cds_end <= tx_end):
        logging.warning(
            "Skipping transcript %s: CDS %d-%d outside transcript bounds %d-%d.",
            gene_id, cds_start, cds_end, tx_start, tx_end,
        )
        return ("parse_error", None)
    if strand not in ("+", "-"):
        logging.warning("Skipping transcript %s: unrecognized strand %r.", gene_id, strand)
        return ("parse_error", None)

    row = extract_utr_bed12(
        gene_id, chrom, strand, cds_start, cds_end, exon_starts, exon_ends, utr_type,
    )
    return ("ok" if row else "no_utr", row)


def _utr_output_path(output_file: str, ut: str) -> str:
    """Derive the per-UTR output path without double suffixes."""
    if output_file.endswith(".bed12"):
        return output_file[: -len(".bed12")] + f"_{ut}.bed12"
    if output_file.endswith(".bed"):
        return output_file[: -len(".bed")] + f"_{ut}.bed"
    return f"{output_file}_{ut}"


def convert(input_file: str, output_file: str, utr_type: str, threads: int = 1) -> None:
    """Convert genePred to BED12 with UTR regions.

    Args:
        input_file: Input genePred file.
        output_file: Output BED12 file (or prefix with --utr both).
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
        out_path = _utr_output_path(output_file, ut) if utr_type == "both" else output_file

        jobs = [(line, ut) for line in lines]
        if threads > 1:
            with Pool(threads) as pool:
                results = pool.map(process_line, jobs)
        else:
            results = [process_line(job) for job in jobs]

        n_written = 0
        n_parse_errors = 0
        n_noncoding = 0
        with open(out_path, "w") as fh:
            for status, row in results:
                if status == "parse_error":
                    n_parse_errors += 1
                elif status == "noncoding":
                    n_noncoding += 1
                elif row:
                    fh.write("\t".join(map(str, row)) + "\n")
                    n_written += 1

        if n_parse_errors:
            logging.warning("%d lines failed to parse for %s.", n_parse_errors, ut)
        if n_noncoding:
            logging.info("Skipped %d non-coding transcripts for %s.", n_noncoding, ut)
        logging.info("%s: %d regions -> %s", ut, n_written, out_path)

        if n_written == 0 and n_parse_errors > 0:
            logging.error("No %s regions extracted (all input lines failed to parse).", ut)
            sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract UTR regions from genePred annotations as BED12 format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python genepred_utr_to_bed12.py -i annotation.gp -o utr.bed12 --utr 5UTR
  python genepred_utr_to_bed12.py -i annotation.gp -o utr.bed12 --utr both -t 4

notes:
  On the minus strand the 5'UTR lies right of cdsEnd and the 3'UTR left
  of cdsStart (plus-strand files keep the conventional sides).
  Non-coding rows (cdsStart == cdsEnd) are skipped. Bin-prefixed
  11-column UCSC knownGene exports are auto-detected.
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
