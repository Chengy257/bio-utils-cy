#!/usr/bin/env python3
#########################################################################
# File Name: bed12_effective_length.py
# Author: ChengYu
# Description: Calculate non-redundant genomic length from BED12 files
#              by merging overlapping regions.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: missing `import os` caused an unconditional NameError whenever
#     pybedtools was importable (previously masked by the dependency
#     check exiting first).
#   - REWRITE: intervals are now merged in pure Python (no pybedtools /
#     bedtools binary needed); input is sorted internally before
#     merging, per chromosome, so unsorted input no longer yields wrong
#     partial merges.
#   - FIX: merged-BED output preserves the remaining BED columns of the
#     first interval in each merged group and is written sorted.
#   - Malformed BED lines are skipped with a warning and counted;
#     exit 1 when no valid interval remains.
#########################################################################
"""Calculate non-redundant genomic length from BED files.

Merges overlapping (and bookended) intervals per chromosome and computes
total effective length. Works with BED6, BED12, and generic BED formats
(only the first three columns are used for merging; other columns are
preserved from the first interval of each merged group).

Requires: nothing beyond the Python standard library.
"""

import argparse
import logging
import sys
from pathlib import Path

__version__ = "1.1.0"


def read_intervals(input_file: str):
    """Read (chrom, start, end, extra_columns) tuples from a BED file.

    Blank lines, comments, and track/browser lines are ignored;
    malformed interval lines are skipped with a warning.
    """
    intervals = []
    n_malformed = 0
    with open(input_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                n_malformed += 1
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                n_malformed += 1
                continue
            if end < start:
                logging.warning("Skipping interval with end < start: %s", line)
                n_malformed += 1
                continue
            intervals.append((fields[0], start, end, fields[3:]))
    if n_malformed:
        logging.warning("Skipped %d malformed BED lines.", n_malformed)
    return intervals


def merge_intervals(intervals):
    """Merge overlapping/bookended intervals per chromosome.

    Args:
        intervals: Iterable of (chrom, start, end, extra) tuples.

    Returns:
        List of (chrom, start, end, extra) with extra taken from the
        first interval of each merged group, sorted by chrom then start.
    """
    by_chrom = {}
    for chrom, start, end, extra in intervals:
        by_chrom.setdefault(chrom, []).append((start, end, extra))

    merged = []
    for chrom in sorted(by_chrom):
        runs = sorted(by_chrom[chrom], key=lambda x: (x[0], x[1]))
        cur_start, cur_end, cur_extra = runs[0]
        for start, end, extra in runs[1:]:
            if start <= cur_end:  # overlapping or bookended
                if end > cur_end:
                    cur_end = end
            else:
                merged.append((chrom, cur_start, cur_end, cur_extra))
                cur_start, cur_end, cur_extra = start, end, extra
        merged.append((chrom, cur_start, cur_end, cur_extra))
    return merged


def calculate_effective_length(input_file: str, output_file: str = None) -> int:
    """Calculate non-redundant genomic length from a BED file.

    Args:
        input_file: Input BED file path.
        output_file: Optional path to write merged BED regions.

    Returns:
        Total non-redundant length.

    Raises:
        ValueError: if the file contains no valid intervals.
    """
    intervals = read_intervals(input_file)
    if not intervals:
        raise ValueError(f"no valid BED intervals found in {input_file}")

    merged = merge_intervals(intervals)
    total_length = sum(end - start for _, start, end, _ in merged)

    if output_file:
        with open(output_file, "w") as fh:
            for chrom, start, end, extra in merged:
                fh.write("\t".join([chrom, str(start), str(end)] + list(extra)) + "\n")
        logging.info("Merged regions saved to %s", output_file)

    return total_length


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate non-redundant genomic length from BED files by merging overlapping regions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Print effective length
  python bed12_effective_length.py -i regions.bed

  # Also save merged regions
  python bed12_effective_length.py -i regions.bed -o merged.bed
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input BED file.")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output merged BED file (optional).")
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

    try:
        length = calculate_effective_length(args.input, args.output)
    except ValueError as e:
        logging.error("%s", e)
        sys.exit(1)

    print(f"Effective length: {length:,} bp")


if __name__ == "__main__":
    main()
