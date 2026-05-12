#!/usr/bin/python
#########################################################################
# File Name: aligned_similarity.py
# Author: Code Copilot
# Description: Analyze aligned FASTA files for similarity metrics directly.
# Created Time: Mon 30 Dec 2024 06:00:00 PM CST
#########################################################################
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

def setup_logger():
    """Configure logger settings."""
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )

def validate_fasta(file_path):
    """
    Validate the input file is in FASTA format and contains valid sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        bool: True if file is valid, False otherwise.
    """
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        for record in records:
            if not record.seq:
                logging.error(f"Invalid sequence in {record.id}")
                return False
        logging.info(f"FASTA format validated for {file_path}.")
        return True
    except Exception as e:
        logging.error(f"Error validating FASTA file: {e}")
        return False

def clean_aligned_sequence(ref_seq, target_seq):
    """
    Remove gaps from aligned sequences, ensuring lengths remain consistent.

    Args:
        ref_seq (Seq): Reference aligned sequence.
        target_seq (Seq): Target aligned sequence.

    Returns:
        tuple: Cleaned reference and target sequences with consistent lengths.
    """
    cleaned_ref = []
    cleaned_target = []

    for ref_base, target_base in zip(ref_seq, target_seq):
        if not (ref_base == "-" and target_base == "-"):
            cleaned_ref.append(ref_base)
            cleaned_target.append(target_base)

    if len(cleaned_ref) != len(cleaned_target):
        logging.warning("Cleaned sequences have inconsistent lengths.")

    return Seq("".join(cleaned_ref)), Seq("".join(cleaned_target))

def calculate_similarity(seq1, seq2):
    """
    Calculate the similarity between two sequences.

    Args:
        seq1 (Seq): First sequence.
        seq2 (Seq): Second sequence.

    Returns:
        dict: Dictionary containing similarity metrics.
    """
    if len(seq1) != len(seq2):
        return {"matches": None, "total": None, "similarity": None, "error": "Aligned sequences have unequal lengths."}

    matches = sum(base1 == base2 for base1, base2 in zip(seq1, seq2))
    total = len(seq1)
    similarity = matches / total if total > 0 else None
    return {"matches": matches, "total": total, "similarity": similarity, "error": None}

def process_aligned_file(file_path, reference_species):
    """
    Process an aligned FASTA file to calculate similarity metrics.

    Args:
        file_path (str): Path to the aligned FASTA file.
        reference_species (str): Reference species name.

    Returns:
        list: List of similarity results for the reference and target species.
    """
    if not validate_fasta(file_path):
        return [{"file": os.path.basename(file_path), "reference": None, "target": None, "matches": None, "total": None, "similarity": None, "error": "Invalid FASTA file."}]

    records = list(SeqIO.parse(file_path, "fasta"))
    ref_record = None

    # Identify reference species
    for record in records:
        if record.id.startswith(reference_species + "."):
            ref_record = record
            break

    if not ref_record:
        return [{"file": os.path.basename(file_path), "reference": None, "target": None, "matches": None, "total": None, "similarity": None, "error": "Reference species not found."}]

    results = []

    # Process each target species
    for record in records:
        if record.id != ref_record.id:
            cleaned_ref, cleaned_target = clean_aligned_sequence(ref_record.seq, record.seq)
            similarity_metrics = calculate_similarity(cleaned_ref, cleaned_target)
            results.append({
                "file": os.path.basename(file_path),
                "reference": ref_record.id,
                "target": record.id,
                **similarity_metrics
            })

    return results

def analyze_aligned_files_concurrent(input_dir, reference_species, output_file, threads):
    """
    Analyze aligned FASTA files for similarity metrics with concurrent processing.

    Args:
        input_dir (str): Input directory with aligned FASTA files.
        reference_species (str): Reference species name.
        output_file (str): Output CSV file for results.
        threads (int): Number of threads for concurrent processing.
    """
    results = []
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".aligned.fa")]

    if not files:
        logging.error("No aligned FASTA files found in the input directory.")
        return

    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_file = {executor.submit(process_aligned_file, file_path, reference_species): file_path for file_path in files}
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                file_results = future.result()
                results.extend(file_results)
            except Exception as e:
                logging.error(f"Error processing file {file_path}: {e}")
                results.append({"file": os.path.basename(file_path), "reference": None, "target": None, "matches": None, "total": None, "similarity": None, "error": str(e)})

    if results:
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        logging.info(f"Results written to {output_file}")
    else:
        logging.warning("No results to write. Check input files and reference species.")

def main():
    setup_logger()
    parser = argparse.ArgumentParser(description="Analyze aligned FASTA files for similarity metrics.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input directory with aligned FASTA files.")
    parser.add_argument("-r", "--reference", type=str, required=True, help="Reference species name.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file for results.")
    parser.add_argument("-p", "--threads", type=int, default=4, help="Number of threads for concurrent processing.")
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        logging.error(f"Input directory not found: {args.input}")
        return

    analyze_aligned_files_concurrent(args.input, args.reference, args.output, args.threads)

if __name__ == "__main__":
    main()
