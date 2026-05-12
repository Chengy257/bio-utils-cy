#!/usr/bin/env python3
#########################################################################
# File Name: reverse_complement.py
# Author: ChengYu
# Description: Compute reverse complement, complement, or reverse of
#              nucleotide sequences from FASTA files or command-line input.
# Created Time: 2026
#########################################################################
"""Compute reverse complement, complement, or reverse of nucleotide sequences.

Supports input from FASTA files or stdin strings. Outputs result to stdout
or a specified file.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict

__version__ = "1.0.0"

COMPLEMENT_MAP: Dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
    "a": "t", "t": "a", "g": "c", "c": "g", "n": "n",
}

VALID_MODES = ("revcomp", "comp", "rev")


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a nucleotide sequence."""
    return "".join(COMPLEMENT_MAP[base] for base in reversed(seq))


def complement(seq: str) -> str:
    """Return complement of a nucleotide sequence."""
    return "".join(COMPLEMENT_MAP[base] for base in seq)


def reverse(seq: str) -> str:
    """Return reversed sequence."""
    return seq[::-1]


MODE_FUNCS = {
    "revcomp": reverse_complement,
    "comp": complement,
    "rev": reverse,
}


def parse_fasta(file_path: str) -> list:
    """Parse a FASTA file into a list of (header, sequence) tuples."""
    records = []
    header = None
    seq_parts = []

    with open(file_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line.split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def validate_sequence(seq: str) -> None:
    """Check that all bases are valid nucleotides."""
    valid = set(COMPLEMENT_MAP.keys()) | {"U", "u", "-", "."}
    invalid = set(seq) - valid
    if invalid:
        logging.warning("Non-standard bases found: %s", invalid)


def process_fasta_file(file_path: str, mode: str, output: str = None) -> None:
    """Process a FASTA file and apply the selected operation."""
    func = MODE_FUNCS[mode]
    records = parse_fasta(file_path)
    logging.info("Read %d sequences from %s", len(records), file_path)

    out_fh = open(output, "w") if output else sys.stdout
    try:
        for header, seq in records:
            validate_sequence(seq)
            result = func(seq)
            out_fh.write(f"{header}\n{result}\n")
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

    logging.info("Done. Mode: %s", mode)


def process_string(seq: str, mode: str) -> None:
    """Process a single sequence string and print the result."""
    validate_sequence(seq)
    func = MODE_FUNCS[mode]
    print(func(seq))


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Compute reverse complement, complement, or reverse of nucleotide sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Reverse complement from FASTA file
  python reverse_complement.py -i genes.fa -m revcomp

  # Complement a sequence string
  python reverse_complement.py -s ATCGATCG -m comp

  # Reverse complement to output file
  python reverse_complement.py -i genes.fa -m revcomp -o output.fa
""",
    )
    parser.add_argument(
        "-i", "--input", type=str, help="Input FASTA file path."
    )
    parser.add_argument(
        "-s", "--seq", type=str, help="Input sequence string (use instead of -i)."
    )
    parser.add_argument(
        "-m", "--mode", type=str, default="revcomp",
        choices=VALID_MODES,
        help="Operation mode (default: revcomp).",
    )
    parser.add_argument(
        "-o", "--output", type=str, default=None,
        help="Output file path (default: stdout).",
    )
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=getattr(logging, args.log_level),
    )

    if not args.input and not args.seq:
        parser.error("Either --input or --seq must be provided.")

    if args.input and args.seq:
        parser.error("Specify only one of --input or --seq, not both.")

    if args.input:
        path = Path(args.input)
        if not path.is_file():
            logging.error("Input file not found: %s", args.input)
            sys.exit(1)
        process_fasta_file(args.input, args.mode, args.output)
    else:
        process_string(args.seq, args.mode)


if __name__ == "__main__":
    main()
