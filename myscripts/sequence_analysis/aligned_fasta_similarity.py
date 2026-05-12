#!/usr/bin/env python3
#########################################################################
# File Name: aligned_fasta_similarity.py
# Author: ChengYu
# Description: Analyze aligned FASTA files for direct pairwise
#              similarity metrics (match count, identity) without BLAST.
# Created Time: 2026
#########################################################################
"""Analyze aligned FASTA files for pairwise similarity without BLAST.

For each multi-species aligned FASTA file, computes match count, total
positions, and similarity percentage between a reference and all other
species by direct character comparison (no external tools needed).
"""

import argparse
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "1.0.0"


def validate_fasta(file_path: str) -> bool:
    """Validate that a file contains parseable FASTA with non-empty sequences."""
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            logging.error("No sequences found in %s", file_path)
            return False
        for rec in records:
            if not rec.seq:
                logging.error("Empty sequence in %s: %s", file_path, rec.id)
                return False
        return True
    except Exception as e:
        logging.error("Error parsing FASTA %s: %s", file_path, e)
        return False


def clean_aligned_pair(ref_seq: Seq, target_seq: Seq) -> tuple:
    """Remove double-gap columns from an aligned pair.

    Positions where both sequences have a gap are excluded.

    Args:
        ref_seq: Reference aligned sequence.
        target_seq: Target aligned sequence.

    Returns:
        Tuple of (cleaned_ref, cleaned_target) as Seq objects.
    """
    if len(ref_seq) != len(target_seq):
        logging.warning("Sequence length mismatch (%d vs %d).", len(ref_seq), len(target_seq))
        min_len = min(len(ref_seq), len(target_seq))
        ref_seq = ref_seq[:min_len]
        target_seq = target_seq[:min_len]

    ref_clean = []
    target_clean = []
    for r, t in zip(str(ref_seq), str(target_seq)):
        if r == "-" and t == "-":
            continue
        ref_clean.append(r)
        target_clean.append(t)

    return Seq("".join(ref_clean)), Seq("".join(target_clean))


def calculate_similarity(seq1: Seq, seq2: Seq) -> dict:
    """Calculate similarity metrics between two equal-length sequences.

    Args:
        seq1: First sequence.
        seq2: Second sequence.

    Returns:
        Dict with 'matches', 'total', 'similarity', 'gaps' keys.
    """
    if len(seq1) != len(seq2):
        logging.warning("Length mismatch in calculate_similarity.")
        return {"matches": None, "total": None, "similarity": None, "gaps": None}

    matches = 0
    gaps = 0
    total = len(seq1)
    for a, b in zip(str(seq1), str(seq2)):
        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            matches += 1

    effective_total = total - gaps
    similarity = matches / effective_total if effective_total > 0 else None

    return {
        "matches": matches,
        "total": total,
        "gaps": gaps,
        "similarity": similarity,
    }


def extract_species(record_id: str) -> str:
    """Extract species name from record ID (text before first dot)."""
    return record_id.split(".")[0]


def process_file(file_path: str, reference_species: str) -> list:
    """Process a single aligned FASTA file.

    Args:
        file_path: Path to aligned FASTA file.
        reference_species: Reference species name.

    Returns:
        List of result dicts.
    """
    if not validate_fasta(file_path):
        return []

    records = list(SeqIO.parse(file_path, "fasta"))

    ref_record = None
    for rec in records:
        if extract_species(rec.id) == reference_species:
            ref_record = rec
            break

    if ref_record is None:
        logging.warning("Reference '%s' not found in %s.", reference_species, file_path)
        return []

    results = []
    for rec in records:
        if rec.id == ref_record.id:
            continue

        clean_ref, clean_target = clean_aligned_pair(ref_record.seq, rec.seq)
        metrics = calculate_similarity(clean_ref, clean_target)

        results.append({
            "file": os.path.basename(file_path),
            "reference": ref_record.id,
            "target": rec.id,
            **metrics,
        })

    return results


def analyze_directory(
    input_dir: str,
    reference_species: str,
    output_file: str,
    threads: int = 4,
    suffix: str = ".aligned.fa",
) -> None:
    """Analyze all aligned FASTA files in a directory.

    Args:
        input_dir: Directory with aligned FASTA files.
        reference_species: Reference species name.
        output_file: Output CSV path.
        threads: Number of parallel workers.
        suffix: File suffix to match.
    """
    files = sorted(
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(suffix)
    )

    if not files:
        logging.error("No files matching '*%s' found in %s", suffix, input_dir)
        sys.exit(1)

    logging.info("Processing %d files with %d threads.", len(files), threads)

    results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_map = {
            executor.submit(process_file, fp, reference_species): fp
            for fp in files
        }
        for future in as_completed(future_map):
            fp = future_map[future]
            try:
                file_results = future.result()
                results.extend(file_results)
                logging.info("Completed: %s (%d pairs)", os.path.basename(fp), len(file_results))
            except Exception as e:
                logging.error("Failed to process %s: %s", fp, e)

    if not results:
        logging.warning("No results produced.")
        return

    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    logging.info("Results (%d rows) written to %s", len(df), output_file)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Analyze aligned FASTA files for direct pairwise similarity (no BLAST needed).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python aligned_fasta_similarity.py -i aligned/ -r human -o similarity.csv
  python aligned_fasta_similarity.py -i aligned/ -r human -o similarity.csv --suffix .fa
""",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input directory containing aligned FASTA files.",
    )
    parser.add_argument(
        "-r", "--reference", type=str, required=True,
        help="Reference species name (matched against record ID before first dot).",
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Output CSV file path.",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4,
        help="Number of parallel threads (default: 4).",
    )
    parser.add_argument(
        "--suffix", type=str, default=".aligned.fa",
        help="File suffix to match in input directory (default: .aligned.fa).",
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

    input_dir = Path(args.input)
    if not input_dir.is_dir():
        logging.error("Input directory not found: %s", args.input)
        sys.exit(1)

    analyze_directory(
        input_dir=str(input_dir),
        reference_species=args.reference,
        output_file=args.output,
        threads=args.threads,
        suffix=args.suffix,
    )


if __name__ == "__main__":
    main()
