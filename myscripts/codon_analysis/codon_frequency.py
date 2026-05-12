#!/usr/bin/env python3
#########################################################################
# File Name: codon_frequency.py
# Author: ChengYu
# Description: Calculate codon usage frequencies from nucleotide FASTA.
# Created Time: 2026
#########################################################################
"""Calculate codon usage frequencies from nucleotide FASTA files.

Parses CDS sequences, counts codon occurrences, and outputs a TSV table
with codon counts and percentage frequencies. Supports sorting and
flexible handling of invalid codons.
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

__version__ = "1.0.0"

VALID_NUCLEOTIDES = {"A", "T", "C", "G", "U"}


def parse_fasta_sequences(file_path: str) -> list:
    """Parse a FASTA file into a list of sequence strings.

    Args:
        file_path: Path to FASTA file.

    Yields:
        Sequence strings.
    """
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


def extract_codons(
    seq: str,
    truncate: bool = False,
    skip_invalid: bool = False,
    uppercase: bool = True,
) -> list:
    """Extract valid codons from a nucleotide sequence.

    Args:
        seq: Nucleotide sequence string.
        truncate: If True, truncate to multiple of 3. If False, error on bad length.
        skip_invalid: If True, skip codons with non-standard bases.
        uppercase: Convert sequence to uppercase.

    Returns:
        List of valid 3-letter codon strings.
    """
    if uppercase:
        seq = seq.upper()

    seq_len = len(seq)
    if seq_len == 0:
        return []

    if seq_len % 3 != 0:
        if truncate:
            seq = seq[: seq_len - (seq_len % 3)]
            logging.debug("Truncated sequence from %d to %d bases.", seq_len, len(seq))
        else:
            logging.error("Sequence length %d is not divisible by 3. Use --truncate to auto-fix.", seq_len)
            sys.exit(1)

    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    valid_codons = []

    for codon in codons:
        if all(n in VALID_NUCLEOTIDES for n in codon):
            valid_codons.append(codon)
        elif not skip_invalid:
            logging.error("Invalid codon '%s' found. Use --skip-invalid to skip.", codon)
            sys.exit(1)
        else:
            logging.debug("Skipped invalid codon: %s", codon)

    return valid_codons


def compute_codon_frequency(
    input_path: str,
    output_path: str,
    sort_mode: str = "alpha",
    truncate: bool = False,
    skip_invalid: bool = False,
    uppercase: bool = True,
) -> None:
    """Compute codon frequencies and write TSV output.

    Args:
        input_path: Input FASTA file.
        output_path: Output TSV file.
        sort_mode: 'alpha' or 'freq'.
        truncate: Truncate non-multiple-of-3 sequences.
        skip_invalid: Skip codons with invalid bases.
        uppercase: Convert to uppercase.
    """
    codon_counts = defaultdict(int)
    total = 0

    for seq in parse_fasta_sequences(input_path):
        codons = extract_codons(seq, truncate=truncate, skip_invalid=skip_invalid, uppercase=uppercase)
        for codon in codons:
            codon_counts[codon] += 1
        total += len(codons)

    if total == 0:
        logging.error("No valid codons found in input file.")
        sys.exit(1)

    frequencies = {codon: (count / total) * 100 for codon, count in codon_counts.items()}

    if sort_mode == "alpha":
        sorted_items = sorted(frequencies.items())
    else:
        sorted_items = sorted(frequencies.items(), key=lambda x: (-x[1], x[0]))

    with open(output_path, "w") as fh:
        fh.write("Codon\tCount\tFrequency(%)\n")
        for codon, freq in sorted_items:
            fh.write(f"{codon}\t{codon_counts[codon]}\t{freq:.4f}\n")

    logging.info("Processed %d codons -> %s", total, output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate codon usage frequencies from nucleotide FASTA sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python codon_frequency.py -i cds.fa -o codon_freq.tsv
  python codon_frequency.py -i cds.fa -o codon_freq.tsv --sort freq --truncate
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input nucleotide FASTA file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument(
        "-s", "--sort", choices=["alpha", "freq"], default="alpha",
        help="Sort order: alpha (alphabetical) or freq (descending frequency).",
    )
    parser.add_argument(
        "--truncate", action="store_true",
        help="Truncate sequences to length divisible by 3.",
    )
    parser.add_argument(
        "--skip-invalid", action="store_true",
        help="Skip codons containing non-standard nucleotides.",
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

    compute_codon_frequency(
        input_path=args.input,
        output_path=args.output,
        sort_mode=args.sort,
        truncate=args.truncate,
        skip_invalid=args.skip_invalid,
    )


if __name__ == "__main__":
    main()
