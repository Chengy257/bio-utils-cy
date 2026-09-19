#!/usr/bin/env python3
#########################################################################
# File Name: calculate_cai.py
# Author: ChengYu
# Description: Calculate Codon Adaptation Index (CAI) for target CDS
#              sequences using a reference set of highly expressed genes.
# Created Time: 2026
#########################################################################
"""Calculate Codon Adaptation Index (CAI) for CDS sequences.

Filters target and reference sequences for valid CDS (length divisible
by 3, standard bases only), computes relative adaptiveness weights from
the reference set, then calculates CAI for each target sequence.

Requires: biopython, CAI
"""

import argparse
import logging
import os
import sys
from pathlib import Path

__version__ = "1.0.0"

STANDARD_BASES = set("ATGC")


def is_valid_cds(seq: str) -> bool:
    """Check if sequence is a valid CDS (multiple of 3, standard bases)."""
    seq = seq.upper()
    if len(seq) == 0 or len(seq) % 3 != 0:
        return False
    return all(b in STANDARD_BASES for b in seq)


def filter_reference(input_fasta: str, output_fasta: str) -> tuple:
    """Filter reference sequences to those with valid CDS.

    Args:
        input_fasta: Input reference FASTA.
        output_fasta: Output filtered FASTA.

    Returns:
        Tuple of (total_count, kept_count).
    """
    from Bio import SeqIO

    total = 0
    kept = 0
    with open(output_fasta, "w") as fh:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            if len(rec.seq) % 3 == 0:
                fh.write(f">{rec.id}\n{str(rec.seq)}\n")
                kept += 1
            else:
                logging.debug("Skipped reference %s: length=%d not divisible by 3.", rec.id, len(rec.seq))

    logging.info("Reference: %d total, %d kept, %d skipped.", total, kept, total - kept)
    return total, kept


def filter_targets(input_fasta: str, output_fasta: str) -> tuple:
    """Filter target sequences for valid CDS.

    Args:
        input_fasta: Input target FASTA.
        output_fasta: Output filtered FASTA.

    Returns:
        Tuple of (total_count, kept_count).
    """
    from Bio import SeqIO

    total = 0
    kept = 0
    with open(output_fasta, "w") as fh:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            seq = str(rec.seq).upper()
            if is_valid_cds(seq):
                fh.write(f">{rec.id}\n{seq}\n")
                kept += 1
            else:
                logging.debug("Skipped target %s: invalid CDS.", rec.id)

    logging.info("Targets: %d total, %d kept, %d skipped.", total, kept, total - kept)
    return total, kept


def calculate_cai(
    input_fasta: str,
    ref_fasta: str,
    output_tsv: str,
    work_dir: str = None,
) -> None:
    """Calculate CAI for target sequences against reference.

    Args:
        input_fasta: Target CDS FASTA.
        ref_fasta: Reference CDS FASTA (highly expressed genes).
        output_tsv: Output TSV with CAI values.
        work_dir: Directory for temporary files (default: same as output).
    """
    from Bio import SeqIO
    from CAI import CAI, relative_adaptiveness

    if work_dir is None:
        work_dir = os.path.dirname(os.path.abspath(output_tsv))

    filtered_ref = os.path.join(work_dir, "_filtered_reference.fasta")
    filtered_input = os.path.join(work_dir, "_filtered_input.fasta")

    ref_total, ref_kept = filter_reference(ref_fasta, filtered_ref)
    tgt_total, tgt_kept = filter_targets(input_fasta, filtered_input)

    if ref_kept == 0:
        logging.error("No valid reference sequences after filtering.")
        sys.exit(1)
    if tgt_kept == 0:
        logging.error("No valid target sequences after filtering.")
        sys.exit(1)

    # Compute relative adaptiveness from reference
    reference = [str(rec.seq).upper() for rec in SeqIO.parse(filtered_ref, "fasta")]
    weights = relative_adaptiveness(sequences=reference)
    logging.info("Computed relative adaptiveness from %d reference sequences.", len(reference))

    # Calculate CAI for each target
    computed = 0
    with open(output_tsv, "w") as fh:
        fh.write("seq_id\tCAI\n")
        for rec in SeqIO.parse(filtered_input, "fasta"):
            seq = str(rec.seq).upper()
            try:
                cai_value = CAI(seq, weights=weights)
                fh.write(f"{rec.id}\t{cai_value:.6f}\n")
                computed += 1
            except Exception as e:
                logging.warning("CAI failed for %s: %s", rec.id, e)

    # Clean up temp files
    for tmp in (filtered_ref, filtered_input):
        if os.path.exists(tmp):
            os.unlink(tmp)

    logging.info("CAI computed for %d/%d target sequences -> %s", computed, tgt_kept, output_tsv)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate Codon Adaptation Index (CAI) for target CDS sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python calculate_cai.py -i target_cds.fa -r ribosomal_genes.fa -o cai_results.tsv
  python calculate_cai.py -i all_cds.fa -r highly_expressed.fa -o cai.tsv
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Target CDS FASTA file.")
    parser.add_argument("-r", "--reference", type=str, required=True, help="Reference CDS FASTA (highly expressed genes).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output TSV file with CAI values.")
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

    for path_label, path_val in [("Input", args.input), ("Reference", args.reference)]:
        if not Path(path_val).is_file():
            logging.error("%s file not found: %s", path_label, path_val)
            sys.exit(1)

    try:
        from CAI import CAI, relative_adaptiveness  # noqa: F401
    except ImportError:
        logging.error("Missing dependency: CAI. Install with: pip install cai")
        sys.exit(1)

    calculate_cai(args.input, args.reference, args.output)


if __name__ == "__main__":
    main()
