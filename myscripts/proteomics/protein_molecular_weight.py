#!/usr/bin/env python3
"""
Calculate protein molecular weights from FASTA.

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, TextIO

from Bio import SeqIO

logger = logging.getLogger(__name__)

# Average isotopic masses of amino-acid residues (monoisotopic, Da)
AA_MASS: Dict[str, float] = {
    "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
    "F": 147.06841, "G": 57.02146,  "H": 137.05891, "I": 113.08406,
    "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
    "P": 97.05276,  "Q": 128.05858, "R": 156.10111, "S": 87.03203,
    "T": 101.04768, "V": 99.06841,  "W": 186.07931, "Y": 163.06333,
}

# Water mass added to convert residue weight to molecular weight
WATER_MASS = 18.01056


def calculate_weight(sequence: str) -> float:
    """Sum residue masses and add water to get the molecular weight.

    Parameters
    ----------
    sequence : str
        Protein sequence (one-letter codes).  Non-standard characters are
        silently skipped.

    Returns
    -------
    float
        Molecular weight in Daltons.
    """
    weight = sum(AA_MASS.get(aa, 0.0) for aa in sequence)
    return weight + WATER_MASS


def process_fasta(
    fasta_path: str,
    output_path: Optional[str] = None,
    delimiter: str = "\t",
) -> None:
    """Read a FASTA file, compute molecular weights, and write results.

    Parameters
    ----------
    fasta_path : str
        Path to the input FASTA file.
    output_path : str or None
        Path to the output TSV file.  If *None*, results are printed to
        stdout.
    delimiter : str
        Column delimiter for the output table.
    """
    fasta = Path(fasta_path)
    if not fasta.is_file():
        logger.error("FASTA file not found: %s", fasta_path)
        sys.exit(1)

    records = list(SeqIO.parse(str(fasta), "fasta"))
    if not records:
        logger.error("No sequences found in %s", fasta_path)
        sys.exit(1)

    logger.info("Loaded %d sequences from %s", len(records), fasta_path)

    header = ["ID", "Weight_Da", "Weight_kDa", "Length", "First_AA", "Sequence"]
    rows: List[Dict[str, str]] = []

    for rec in records:
        seq = str(rec.seq).upper()
        valid_len = sum(1 for aa in seq if aa in AA_MASS)
        mw = calculate_weight(seq)
        rows.append({
            "ID": rec.id,
            "Weight_Da": f"{mw:.2f}",
            "Weight_kDa": f"{mw / 1000:.2f}",
            "Length": str(valid_len),
            "First_AA": seq[0] if seq else "",
            "Sequence": seq.strip("."),
        })

    if output_path:
        with open(output_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=header, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
        logger.info("Results written to: %s (%d entries)", output_path, len(rows))
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=header, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="protein_molecular_weight.py",
        description="Calculate protein molecular weights from a FASTA file.",
        epilog=(
            "Examples:\n"
            "  %(prog)s -i proteins.fasta\n"
            "  %(prog)s -i proteins.fasta -o weights.tsv\n"
            "  %(prog)s -i input.fa -o output.tsv --log-level DEBUG\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("-o", "--output", default=None, help="Output TSV file (default: stdout).")
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    process_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
