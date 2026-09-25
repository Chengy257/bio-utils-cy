#!/usr/bin/env python3
#########################################################################
# File Name: amino_acid_frequency.py
# Author: ChengYu
# Description: Calculate amino acid composition frequencies from
#              protein FASTA files.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-25
#   - FIX: with --standard-only the frequencies are now relative to the
#     standard residues only (they sum to 100%); previously non-standard
#     characters (X, *, ...) were excluded from the counts but still
#     inflated the denominator, so the table silently summed to < 100%.
#     The number of skipped residues is logged.
#   - CLEANUP: removed the dead internal `uppercase` parameter (sequences
#     are always uppercased before counting, unchanged behavior).
#########################################################################
"""Calculate amino acid composition from protein FASTA files.

Counts amino acid occurrences and outputs a TSV with counts and
percentage frequencies. Supports sorting and case normalization.
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

__version__ = "1.1.0"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta_sequences(file_path: str):
    """Parse a FASTA file, yielding sequence strings."""
    with open(file_path, "r") as fh:
        seq_parts = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    yield "".join(seq_parts)
                seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            yield "".join(seq_parts)


def compute_aa_frequency(
    input_path: str,
    output_path: str,
    sort_mode: str = "alpha",
    standard_only: bool = False,
) -> None:
    """Compute amino acid frequencies and write TSV output.

    Args:
        input_path: Input protein FASTA file.
        output_path: Output TSV file.
        sort_mode: 'alpha' or 'freq'.
        standard_only: Only count the 20 standard amino acids; frequencies
            are relative to standard residues only.
    """
    counts = defaultdict(int)
    total = 0
    skipped = 0

    for seq in parse_fasta_sequences(input_path):
        seq = seq.upper()
        for aa in seq:
            if standard_only and aa not in STANDARD_AA:
                skipped += 1
                continue
            counts[aa] += 1
            total += 1

    if total == 0:
        logging.error("No amino acids found in input file.")
        sys.exit(1)

    if skipped:
        logging.info("Skipped %d non-standard residue(s)%s.", skipped,
                     " (--standard-only)" if standard_only else "")

    frequencies = {aa: (count / total) * 100 for aa, count in counts.items()}

    if sort_mode == "alpha":
        sorted_items = sorted(frequencies.items())
    else:
        sorted_items = sorted(frequencies.items(), key=lambda x: (-x[1], x[0]))

    with open(output_path, "w") as fh:
        fh.write("Amino_acid\tCount\tFrequency(%)\n")
        for aa, freq in sorted_items:
            fh.write(f"{aa}\t{counts[aa]}\t{freq:.4f}\n")

    logging.info("Processed %d amino acids -> %s", total, output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate amino acid composition frequencies from protein FASTA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python amino_acid_frequency.py -i proteins.fa -o aa_freq.tsv
  python amino_acid_frequency.py -i proteins.fa -o aa_freq.tsv --sort freq --standard-only

notes:
  Sequences are uppercased before counting. Without --standard-only every
  residue character (including *, X, B, Z...) is counted as-is; with
  --standard-only non-standard residues are dropped and frequencies are
  relative to the standard residues only (sum to 100%).
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input protein FASTA file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument(
        "-s", "--sort", choices=["alpha", "freq"], default="alpha",
        help="Sort order: alpha or freq (default: alpha).",
    )
    parser.add_argument(
        "--standard-only", action="store_true",
        help="Count only the 20 standard amino acids (ACDEFGHIKLMNPQRSTVWY); "
             "frequencies are relative to standard residues only.",
    )
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

    compute_aa_frequency(
        input_path=args.input,
        output_path=args.output,
        sort_mode=args.sort,
        standard_only=args.standard_only,
    )


if __name__ == "__main__":
    main()
