#!/usr/bin/python
#########################################################################
# File Name: align_analysis.py
# Author: Code Copilot
# Description: Analyze DNA files using BLAST for similarity metrics.
# Created Time: Sun 29 Dec 2024 04:27:50 PM CST
#########################################################################
import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from tempfile import NamedTemporaryFile
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from Bio.Seq import Seq


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
            if not record.seq or not str(record.seq).isalpha():
                logging.error(f"Invalid sequence in {record.id}")
                return False
        logging.info(f"FASTA format validated for {file_path}.")
        return True
    except Exception as e:
        logging.error(f"Error validating FASTA file: {e}")
        return False


def clean_sequence(seq):
    """
    Remove gaps and ensure sequence length is suitable for translation.

    Args:
        seq (Seq): Input sequence.

    Returns:
        Seq or None: Cleaned sequence or None if invalid.
    """
    cleaned_seq = str(seq).replace("-", "").upper()
    if len(cleaned_seq) < 3:
        logging.warning("Sequence is too short after cleaning. Skipping.")
        return None
    if len(cleaned_seq) % 3 != 0:
        cleaned_seq = cleaned_seq[:-(len(cleaned_seq) % 3)]
    return Seq(cleaned_seq) if cleaned_seq else None


def run_blast(query_seq, target_seq, is_protein=False):
    """
    Run BLAST locally between two sequences with parameters optimized for short sequences.

    Args:
        query_seq (Seq): Query sequence.
        target_seq (Seq): Target sequence.
        is_protein (bool): Whether the sequences are protein.

    Returns:
        dict: Dictionary with BLAST results (e-value, bit score, identity).
    """
    if len(query_seq) < 3 or len(target_seq) < 3:
        logging.warning("Query or target sequence is too short for BLAST. Skipping.")
        return {"e-value": None, "bit score": None, "identity": None}

    with NamedTemporaryFile(delete=False, suffix=".fa") as query_file, NamedTemporaryFile(delete=False, suffix=".fa") as target_file:
        query_file.write(f">query\n{query_seq}\n".encode())
        target_file.write(f">target\n{target_seq}\n".encode())
        query_file.close()
        target_file.close()

        blast_program = "blastp" if is_protein else "blastn"
        result_file_path = NamedTemporaryFile(delete=False, suffix=".xml").name

        blast_options = [
            "-query", query_file.name,
            "-subject", target_file.name,
            "-outfmt", "5",
            "-out", result_file_path
        ]

        if is_protein:
            blast_options.extend(["-matrix", "PAM30", "-seg", "no", "-gapopen", "9", "-gapextend", "1"])
        else:
            blast_options.extend(["-task", "blastn-short", "-dust", "no", "-gapopen", "4", "-gapextend", "2"])

        try:
            subprocess.run(
                [blast_program] + blast_options,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            with open(result_file_path) as blast_out:
                blast_record = NCBIXML.read(blast_out)

            if blast_record.alignments:
                alignment = blast_record.alignments[0]
                hsp = alignment.hsps[0]
                return {"e-value": hsp.expect, "bit score": hsp.bits, "identity": hsp.identities / hsp.align_length}
            else:
                logging.warning("No BLAST alignment found.")
                return {"e-value": None, "bit score": None, "identity": None}
        except subprocess.CalledProcessError as e:
            logging.error(f"BLAST execution failed. Error:\n{e.stderr.decode()}")
            return {"e-value": None, "bit score": None, "identity": None}
        finally:
            os.unlink(query_file.name)
            os.unlink(target_file.name)
            os.unlink(result_file_path)


def process_file(file_path, reference_species):
    """
    Process a single DNA FASTA file, performing BLAST-based similarity calculations.

    Args:
        file_path (str): Path to the DNA FASTA file.
        reference_species (str): Reference species name.

    Returns:
        list: List of dictionaries with similarity results for each species.
    """
    if not validate_fasta(file_path):
        return []

    results = []
    alignment = list(SeqIO.parse(file_path, "fasta"))

    # 提取参考序列
    ref_record = None
    for record in alignment:
        species_name = record.id.split('.')[0]  # 提取物种名称
        if species_name == reference_species:
            ref_record = record
            break

    if not ref_record:
        logging.warning(f"Reference species {reference_species} not found in {file_path}. Skipping.")
        return results

    # 清理参考序列
    ref_record.seq = clean_sequence(ref_record.seq)
    if ref_record.seq is None:
        return results

    # 对比每个物种的序列
    for record in alignment:
        species_name = record.id.split('.')[0]  # 提取物种名称
        if species_name != reference_species:
            record.seq = clean_sequence(record.seq)
            if record.seq is None:
                continue

            # 运行 BLAST
            print(ref_record.seq)
            print(record.seq)
            nucleotide_blast = run_blast(ref_record.seq, record.seq, is_protein=False)
            ref_protein = ref_record.seq.translate()
            target_protein = record.seq.translate()
            protein_blast = run_blast(ref_protein, target_protein, is_protein=True)

            results.append({
                "File": os.path.basename(file_path),
                "Reference": ref_record.id,
                "Species": record.id,
                "Nucleotide e-value": nucleotide_blast["e-value"],
                "Nucleotide bit score": nucleotide_blast["bit score"],
                "Nucleotide identity": nucleotide_blast["identity"],
                "Protein identity": protein_blast["identity"]
            })
    return results



def analyze_dna_files_concurrent(dna_dir, reference_species, output_file, threads):
    """
    Analyze DNA files using BLAST for similarity metrics with concurrent processing.

    Args:
        dna_dir (str): Input directory with DNA FASTA files.
        reference_species (str): Reference species name.
        output_file (str): Output CSV file for results.
        threads (int): Number of threads for concurrent processing.
    """
    results = []
    files = [os.path.join(dna_dir, f) for f in os.listdir(dna_dir) if f.endswith(".dna.fa")]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_file = {executor.submit(process_file, file_path, reference_species): file_path for file_path in files}
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                results.extend(future.result())
            except Exception as e:
                logging.error(f"Error processing file {file_path}: {e}")

    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    logging.info(f"Results written to {output_file}")


def main():
    setup_logger()
    parser = argparse.ArgumentParser(description="Analyze DNA files using BLAST for similarity metrics.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input directory with DNA FASTA files.")
    parser.add_argument("-r", "--reference", type=str, required=True, help="Reference species name.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file for results.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for concurrent processing.")
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        logging.error(f"Input directory not found: {args.input}")
        return

    analyze_dna_files_concurrent(args.input, args.reference, args.output, args.threads)


if __name__ == "__main__":
    main()
