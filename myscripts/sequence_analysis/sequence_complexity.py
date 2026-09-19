#!/usr/bin/env python3
#########################################################################
# File Name: sequence_complexity.py
# Author: ChengYu
# Description: Calculate sequence complexity metrics (block entropy and
#              conditional entropy) for FASTA sequences using sliding
#              window k-mer analysis.
# Created Time: 2026
#########################################################################
"""Calculate sequence complexity metrics for FASTA sequences.

Computes block entropy and conditional entropy using k-mer analysis,
providing measures of sequence randomness and predictability.

Requires: seq_entropies package
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Tuple

__version__ = "1.0.0"


def calc_complexity(seq: str, window: int = 3, base: int = 2) -> Tuple[float, float]:
    """Calculate block entropy and conditional entropy for a sequence.

    Args:
        seq: Input nucleotide or amino acid sequence.
        window: K-mer window size (default: 3).
        base: Logarithm base for entropy calculation (default: 2 for bits).

    Returns:
        Tuple of (block_entropy, conditional_entropy).
    """
    import seq_entropies as se

    symbols = list(seq)
    block_ent = se.block_entropy(symbols, window=window, method="MLE", base=base)
    cond_ent = se.block_cond_entropy(symbols, window=window, method="MLE", base=base)
    return block_ent, cond_ent


def process_fasta(
    input_path: str,
    output_path: str,
    window: int = 3,
    base: int = 2,
) -> None:
    """Process a FASTA file and compute complexity for each sequence.

    Args:
        input_path: Path to input FASTA file.
        output_path: Path to output TSV file.
        window: K-mer window size.
        base: Logarithm base.
    """
    from Bio import SeqIO

    results = []
    total = 0
    for rec in SeqIO.parse(input_path, "fasta"):
        seq = str(rec.seq)
        if not seq:
            logging.warning("Empty sequence: %s, skipping.", rec.id)
            continue
        h, c = calc_complexity(seq, window=window, base=base)
        results.append((rec.id, h, c))
        total += 1

    with open(output_path, "w") as fout:
        fout.write("seq_id\tBlockEntropy\tCondEntropy\n")
        for seq_id, h, c in results:
            fout.write(f"{seq_id}\t{h:.6f}\t{c:.6f}\n")

    logging.info("Processed %d sequences -> %s", total, output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate sequence complexity metrics (block entropy and conditional entropy).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # Basic usage
  python sequence_complexity.py -i sequences.fa -o complexity.tsv

  # Custom window size and base
  python sequence_complexity.py -i sequences.fa -o complexity.tsv -w 4 -b 2
""",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input FASTA file.",
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Output TSV file with complexity metrics.",
    )
    parser.add_argument(
        "-w", "--window", type=int, default=3,
        help="K-mer window size for entropy calculation (default: 3).",
    )
    parser.add_argument(
        "-b", "--base", type=int, default=2,
        help="Logarithm base: 2 for bits, 10 for dits, e for nats (default: 2).",
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

    input_path = Path(args.input)
    if not input_path.is_file():
        logging.error("Input file not found: %s", args.input)
        sys.exit(1)

    if args.window < 1:
        logging.error("Window size must be >= 1.")
        sys.exit(1)

    try:
        import seq_entropies  # noqa: F401
    except ImportError:
        logging.error("Missing dependency: seq_entropies. Install with: pip install seq-entropies")
        sys.exit(1)

    process_fasta(args.input, args.output, window=args.window, base=args.base)


if __name__ == "__main__":
    main()
