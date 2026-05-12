#!/usr/bin/env python3
#########################################################################
# File Name: bed12_effective_length.py
# Author: ChengYu
# Description: Calculate non-redundant genomic length from BED12 files
#              by merging overlapping regions.
# Created Time: 2026
#########################################################################
"""Calculate non-redundant genomic length from BED files.

Merges overlapping intervals and computes total effective length.
Supports BED6, BED12, and generic BED formats.

Requires: pybedtools (or bedtools in PATH)
"""

import argparse
import logging
import sys
from pathlib import Path

__version__ = "1.0.0"


def calculate_effective_length(input_file: str, output_file: str = None) -> int:
    """Calculate non-redundant genomic length from a BED file.

    Args:
        input_file: Input BED file path.
        output_file: Optional path to write merged BED regions.

    Returns:
        Total non-redundant length.
    """
    import pybedtools

    bed = pybedtools.BedTool(input_file)
    merged = bed.merge()

    total_length = 0
    for interval in merged:
        total_length += interval.end - interval.start

    if output_file:
        merged.saveas(output_file)
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
        import pybedtools  # noqa: F401
    except ImportError:
        logging.error("Missing dependency: pybedtools. Install with: pip install pybedtools")
        sys.exit(1)

    length = calculate_effective_length(args.input, args.output)
    print(f"Effective length: {length:,} bp")
    logging.info("Effective length: %d bp", length)


if __name__ == "__main__":
    main()
