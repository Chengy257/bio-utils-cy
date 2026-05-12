#!/usr/bin/env python3
#########################################################################
# File Name: calculate_dnds.py
# Author: ChengYu
# Description: Calculate dN/dS ratios for all species pairs from
#              multi-species aligned DNA FASTA files.
# Created Time: 2026
#########################################################################
"""Calculate dN/dS ratios for all species pairs from aligned DNA FASTA.

Processes directory of .dna.fa files, each containing aligned CDS from
multiple species. Computes dN (non-synonymous) and dS (synonymous)
substitutions for every species pair and outputs a summary TSV.

Requires: biopython
"""

import argparse
import logging
import os
import sys
from itertools import combinations
from pathlib import Path

__version__ = "1.0.0"


def calculate_dn_ds(dna_seq1: str, dna_seq2: str) -> tuple:
    """Calculate dN and dS between two aligned DNA sequences.

    Counts codon positions where the two sequences differ:
    - dN: codon change produces a different amino acid (non-synonymous)
    - dS: codon change produces the same amino acid (synonymous)

    Args:
        dna_seq1: First DNA sequence string.
        dna_seq2: Second DNA sequence string.

    Returns:
        Tuple of (dN, dS, dN/dS ratio).
    """
    from Bio.Data.CodonTable import standard_dna_table

    if len(dna_seq1) != len(dna_seq2):
        raise ValueError("Sequences have different lengths.")

    codons1 = [dna_seq1[i:i + 3] for i in range(0, len(dna_seq1), 3)]
    codons2 = [dna_seq2[i:i + 3] for i in range(0, len(dna_seq2), 3)]

    dn = 0
    ds = 0

    for codon1, codon2 in zip(codons1, codons2):
        if codon1 == codon2:
            continue
        if "-" in codon1 or "-" in codon2 or len(codon1) != 3 or len(codon2) != 3:
            continue

        aa1 = standard_dna_table.forward_table.get(codon1)
        aa2 = standard_dna_table.forward_table.get(codon2)
        if aa1 is None or aa2 is None:
            continue

        if aa1 == aa2:
            ds += 1
        else:
            dn += 1

    ratio = dn / ds if ds > 0 else float("inf")
    return dn, ds, ratio


def process_file(dna_file: str) -> list:
    """Process a single multi-species DNA FASTA file.

    Args:
        dna_file: Path to DNA FASTA file with aligned sequences.

    Returns:
        List of result tuples (sp1, sp2, dN, dS, dN/dS).
    """
    from Bio.SeqIO import parse as seqio_parse

    records = {rec.id: str(rec.seq).upper() for rec in seqio_parse(dna_file, "fasta")}
    results = []

    for sp1, sp2 in combinations(records.keys(), 2):
        seq1 = records[sp1]
        seq2 = records[sp2]

        # Validate: translate and check length consistency
        from Bio.Seq import Seq
        try:
            prot1 = Seq(seq1.replace("-", "")).translate(to_stop=True)
            prot2 = Seq(seq2.replace("-", "")).translate(to_stop=True)
            if len(prot1) != len(prot2):
                logging.debug(
                    "Protein length mismatch for %s vs %s (%d vs %d), skipping.",
                    sp1, sp2, len(prot1), len(prot2),
                )
                continue
        except Exception as e:
            logging.warning("Translation error for %s/%s: %s", sp1, sp2, e)
            continue

        try:
            dn, ds, ratio = calculate_dn_ds(seq1, seq2)
            results.append((sp1, sp2, dn, ds, ratio))
        except ValueError as e:
            logging.warning("dN/dS failed for %s/%s: %s", sp1, sp2, e)

    return results


def batch_process(input_dir: str, output_file: str, suffix: str = ".dna.fa") -> None:
    """Batch process all DNA FASTA files in a directory.

    Args:
        input_dir: Directory with .dna.fa files.
        output_file: Output TSV path.
        suffix: File suffix to match.
    """
    files = sorted(f for f in os.listdir(input_dir) if f.endswith(suffix))
    if not files:
        logging.error("No files matching '*%s' in %s", suffix, input_dir)
        sys.exit(1)

    logging.info("Processing %d files.", len(files))

    with open(output_file, "w") as fh:
        fh.write("Gene\tSpecies1\tSpecies2\tdN\tdS\tdN_dS\n")

        total_pairs = 0
        for fname in files:
            fpath = os.path.join(input_dir, fname)
            gene_name = os.path.splitext(fname)[0]
            logging.info("Processing: %s", fname)

            results = process_file(fpath)
            for sp1, sp2, dn, ds, ratio in results:
                ratio_str = f"{ratio:.4f}" if ratio != float("inf") else "Inf"
                fh.write(f"{gene_name}\t{sp1}\t{sp2}\t{dn}\t{ds}\t{ratio_str}\n")
            total_pairs += len(results)

    logging.info("Completed: %d pairs -> %s", total_pairs, output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate dN/dS ratios for all species pairs from aligned DNA FASTA files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv
  python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv --suffix .fa
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input directory with multi-species DNA FASTA files.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file.")
    parser.add_argument("--suffix", type=str, default=".dna.fa", help="File suffix to match (default: .dna.fa).")
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

    if not Path(args.input).is_dir():
        logging.error("Input directory not found: %s", args.input)
        sys.exit(1)

    batch_process(args.input, args.output, args.suffix)


if __name__ == "__main__":
    main()
