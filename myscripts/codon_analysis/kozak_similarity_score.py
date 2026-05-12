#!/usr/bin/env python3
#########################################################################
# File Name: kozak_similarity_score.py
# Author: ChengYu
# Description: Calculate Kozak consensus sequence similarity score
#              for translation initiation sites.
# Created Time: 2026
#########################################################################
"""Calculate Kozak consensus sequence similarity score.

Scores 23-bp sequences centered on the start codon (ATG at positions
10-12) using a position-specific scoring matrix derived from annotated
TIS data. Supports single sequence input or batch FASTA processing.

Reference: https://github.com/Agleason1/TIS-Predictor
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

__version__ = "1.0.0"

# Position-specific scoring matrix (23 positions x 5 bases: A, T, G, C, N)
# Columns: A=0, T=1, G=2, C=3, N=4
KOZAK_WEIGHTS = np.array([
    [0.04210526, 0.00000000, 0.03157895, 0.05263158, 0.00000000],
    [0.04210526, 0.05263158, 0.10526316, 0.06250000, 0.00000000],
    [0.03157895, 0.04210526, 0.05263158, 0.07368421, 0.00000000],
    [0.03157895, 0.01052632, 0.04210526, 0.05263158, 0.00000000],
    [0.08421053, 0.07368421, 0.18947368, 0.10526316, 0.00000000],
    [0.04210526, 0.05263158, 0.05263158, 0.08421053, 0.00000000],
    [0.12631579, 0.06250000, 0.12631579, 0.21052632, 0.00000000],
    [0.83157895, 0.12631579, 0.65263158, 0.16842105, 0.00000000],
    [0.15789474, 0.06315789, 0.11578947, 0.20000000, 0.00000000],
    [0.21052632, 0.09473684, 0.31578947, 0.51578947, 0.00000000],
    [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
    [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
    [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
    [0.24210526, 0.16666667, 0.53684211, 0.13684211, 0.00000000],
    [0.15789474, 0.09473684, 0.09473684, 0.24210526, 0.00000000],
    [0.05263158, 0.08421053, 0.14736842, 0.09473684, 0.00000000],
    [0.07216495, 0.05263158, 0.10526316, 0.06315789, 0.00000000],
    [0.00000000, 0.00000000, 0.00000000, 0.05263158, 0.00000000],
    [0.05263158, 0.05263158, 0.10526316, 0.09473684, 0.00000000],
    [0.04210526, 0.03157895, 0.05263158, 0.04210526, 0.00000000],
    [0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000],
    [0.04210526, 0.04210526, 0.08421053, 0.07368421, 0.00000000],
    [0.06250000, 0.04210526, 0.09473684, 0.05263158, 0.00000000],
])

BASE_TO_IDX = {"A": 0, "T": 1, "G": 2, "C": 3}
REQUIRED_LENGTH = 23


def kozak_score(sequence: str) -> float:
    """Calculate normalized Kozak similarity score for a 23-bp sequence.

    The ATG start codon should be centered at positions 10-12 (0-indexed).
    Score is normalized by the maximum possible score (range 0-1).

    Args:
        sequence: 23-bp nucleotide string centered on ATG.

    Returns:
        Normalized score between 0 and 1.

    Raises:
        ValueError: If sequence is not 23 bases long.
    """
    if len(sequence) != REQUIRED_LENGTH:
        raise ValueError(f"Sequence must be {REQUIRED_LENGTH} bases long, got {len(sequence)}.")

    seq = sequence.upper().replace("U", "T")
    score = 0.0
    for i, base in enumerate(seq):
        idx = BASE_TO_IDX.get(base, 4)  # Unknown bases map to N column
        score += KOZAK_WEIGHTS[i][idx]

    max_score = KOZAK_WEIGHTS.max(axis=1).sum()
    return score / max_score


def process_fasta(input_path: str, output_path: str) -> None:
    """Process a FASTA file and score each sequence.

    Args:
        input_path: Input FASTA file with 23-bp sequences.
        output_path: Output TSV file.
    """
    results = []
    header = None
    seq_parts = []

    with open(input_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts)
                    try:
                        score = kozak_score(seq)
                        results.append((header, f"{score:.6f}"))
                    except ValueError as e:
                        logging.warning("Skipping %s: %s", header, e)
                header = line.split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    # Last record
    if header is not None:
        seq = "".join(seq_parts)
        try:
            score = kozak_score(seq)
            results.append((header, f"{score:.6f}"))
        except ValueError as e:
            logging.warning("Skipping %s: %s", header, e)

    with open(output_path, "w") as fh:
        fh.write("seq_id\tKozak_score\n")
        for seq_id, score in results:
            fh.write(f"{seq_id}\t{score}\n")

    logging.info("Scored %d sequences -> %s", len(results), output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate Kozak consensus sequence similarity score for TIS sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
The input sequence must be exactly 23 bp with the ATG start codon centered
at positions 10-12 (0-indexed). The score ranges from 0 to 1.

examples:
  # Score a single sequence
  python kozak_similarity_score.py -s GCCACCATGGCGATCGATCGATC

  # Batch score from FASTA
  python kozak_similarity_score.py -i tis_sequences.fa -o scores.tsv
""",
    )
    parser.add_argument("-i", "--input", type=str, help="Input FASTA file with 23-bp TIS sequences.")
    parser.add_argument("-o", "--output", type=str, help="Output TSV file (required with -i).")
    parser.add_argument("-s", "--seq", type=str, help="Single 23-bp sequence to score.")
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

    if not args.input and not args.seq:
        parser.error("Either --input or --seq must be provided.")
    if args.input and args.seq:
        parser.error("Specify only one of --input or --seq.")

    if args.seq:
        try:
            score = kozak_score(args.seq)
            print(f"Kozak similarity score: {score:.6f}")
        except ValueError as e:
            logging.error(str(e))
            sys.exit(1)
    else:
        if not args.output:
            parser.error("--output is required when using --input.")
        if not Path(args.input).is_file():
            logging.error("Input file not found: %s", args.input)
            sys.exit(1)
        process_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
