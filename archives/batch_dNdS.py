#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author:       ChengYu
# Description:  Calculate dN/dS for all species pairs using dna.fa files only.
# Created:      2022/12/07
# Last Edited:  2024/11/21

import os
import argparse
from itertools import combinations
from Bio.Seq import Seq
from Bio.SeqIO import parse
from Bio.Data.CodonTable import standard_dna_table


def calculate_dn_ds(dna_seq1, dna_seq2):
    """
    Calculate dN/dS (non-synonymous/synonymous substitution ratio) between two DNA sequences.

    Args:
        dna_seq1 (str): First DNA sequence.
        dna_seq2 (str): Second DNA sequence.

    Returns:
        tuple: dN (non-synonymous substitutions), dS (synonymous substitutions), dN/dS ratio.
    """
    codons1 = [str(dna_seq1[i:i+3]) for i in range(0, len(dna_seq1), 3)]
    codons2 = [str(dna_seq2[i:i+3]) for i in range(0, len(dna_seq2), 3)]

    if len(codons1) != len(codons2):
        raise ValueError("Sequences have different lengths, cannot calculate dN/dS.")

    dN = 0
    dS = 0

    for codon1, codon2 in zip(codons1, codons2):
        if codon1 == codon2:
            continue  # No substitution
        if "-" in codon1 or "-" in codon2 or len(codon1) != 3 or len(codon2) != 3:
            continue  # Skip invalid codons
        aa1 = standard_dna_table.forward_table.get(codon1, None)
        aa2 = standard_dna_table.forward_table.get(codon2, None)
        if aa1 is None or aa2 is None:
            continue  # Skip invalid codons
        if aa1 == aa2:
            dS += 1  # Synonymous substitution
        else:
            dN += 1  # Non-synonymous substitution

    dn_ds_ratio = dN / dS if dS > 0 else float('inf')
    return dN, dS, dn_ds_ratio


def process_dna_file(dna_file):
    """
    Calculate dN/dS for all species pairs in a single dna.fa file.

    Args:
        dna_file (str): Path to the DNA FASTA file.

    Returns:
        list: List of tuples containing results for each species pair.
    """
    results = []
    dna_records = {rec.id: rec.seq for rec in parse(dna_file, "fasta")}
    species_pairs = combinations(dna_records.keys(), 2)

    for species1, species2 in species_pairs:
        dna_seq1 = dna_records.get(species1)
        dna_seq2 = dna_records.get(species2)

        if not dna_seq1 or not dna_seq2:
            print(f"Warning: Missing sequences for {species1} or {species2}, skipping.")
            continue

        # Ensure DNA translates correctly to protein
        try:
            protein_seq1 = dna_seq1.replace("-", "").translate(to_stop=True)
            protein_seq2 = dna_seq2.replace("-", "").translate(to_stop=True)

            # Check for mismatches in translation length
            if len(protein_seq1) != len(protein_seq2):
                print(f"Warning: Protein mismatch for {species1} and {species2}, skipping.")
                continue
        except Exception as e:
            print(f"Error translating sequences for {species1} or {species2}: {e}")
            continue

        # Calculate dN/dS
        try:
            dN, dS, dn_ds_ratio = calculate_dn_ds(dna_seq1, dna_seq2)
            results.append((species1, species2, dN, dS, dn_ds_ratio))
        except ValueError as e:
            print(f"Error calculating dN/dS for {species1} and {species2}: {e}")
    return results


def batch_process(input_dir, output_file):
    """
    Batch process all dna.fa files in the input directory to calculate dN/dS.

    Args:
        input_dir (str): Path to the directory containing dna.fa files.
        output_file (str): Path to the output TSV file.

    Returns:
        None
    """
    with open(output_file, "w") as out:
        out.write("File_Group\tSpecies1\tSpecies2\tdN\tdS\tdN/dS\n")

        dna_files = [f for f in os.listdir(input_dir) if f.endswith(".dna.fa")]

        for dna_file in dna_files:
            dna_path = os.path.join(input_dir, dna_file)
            print(f"Processing: {dna_file}...")

            results = process_dna_file(dna_path)

            group = os.path.splitext(dna_file)[0]
            for species1, species2, dN, dS, dn_ds_ratio in results:
                out.write(f"{group}\t{species1}\t{species2}\t{dN}\t{dS}\t{dn_ds_ratio:.4f}\n")


def main():
    parser = argparse.ArgumentParser(description="Batch calculate dN/dS for all dna.fa files.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input directory containing dna.fa files.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output TSV file.")

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"Error: Input directory not found: {args.input}")
        return

    try:
        batch_process(args.input, args.output)
        print(f"Completed. dN/dS results saved to {args.output}")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
