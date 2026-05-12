#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author:       ChengYu
# Descriptions: Extract subregions from multi-alignment MAF files and convert to FASTA.
#               Supports ungapped and translated sequences. MAF files must be split by chromosome.
# Created:      2022/12/07
# Last Edited:  2024/11/21

import os
import argparse
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import SeqIO

def maf_parse(maf_dir, gene_pred, species):
    """
    Extract and process regions from a MAF file based on genePred data.

    Args:
        maf_dir (str): Directory containing chromosome-split MAF files.
        gene_pred (list): Parsed gene annotation data for a single entry.
        species (str): Species name used in the MAF file.
    """
    try:
        # Extract required fields from gene_pred
        name, chr_name, strand, cds_start, cds_end = gene_pred[0], gene_pred[1], gene_pred[2], gene_pred[8], gene_pred[9]
        strand = 1 if strand == "+" else -1

        # Prepare paths and filenames
        maf_file = os.path.join(maf_dir, f"{chr_name}.maf")
        maf_index = f"{maf_file}.index"
        output_dna = f"{name}.dna.fa"
        output_protein = f"{name}.pep.fa"
        ref = f"{species}.{chr_name}"

        # Check if MAF file exists
        if not os.path.isfile(maf_file):
            raise FileNotFoundError(f"MAF file not found: {maf_file}")

        # Parse start and end positions
        start_positions = list(map(int, cds_start.strip(",").split(",")))
        end_positions = list(map(int, cds_end.strip(",").split(",")))

        # Index and extract spliced regions
        idx = MafIO.MafIndex(maf_index, maf_file, ref)
        region = idx.get_spliced(start_positions, end_positions, strand=strand)
        
        # Write extracted region to DNA FASTA
        AlignIO.write(region, output_dna, "fasta")

        # Process DNA sequences: ungap and translate to protein
        with open(output_dna, "r") as dna_in, open(output_protein, "w") as protein_out:
            for seq_record in SeqIO.parse(dna_in, "fasta"):
                protein_out.write(f">{seq_record.id}\n")
                seq_ungapped = seq_record.seq.replace("-","")
                protein_out.write(f"{seq_ungapped.translate()}\n")
    except Exception as e:
        print(f"Error processing gene_pred entry {gene_pred}: {e}")

def parse_gene_pred(file_path):
    """
    Parse the genePred file into a list of entries.

    Args:
        file_path (str): Path to the genePred file.

    Returns:
        list: Parsed genePred data entries.
    """
    try:
        with open(file_path, "r") as f:
            return [line.strip().split("\t") for line in f if line.strip()]
    except Exception as e:
        raise IOError(f"Error reading genePred file {file_path}: {e}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract subregions from MAF files and output as FASTA.")
    parser.add_argument("-g", "--genePred", type=str, required=True, help="Path to the gene annotation file (genePred format).")
    parser.add_argument("-d", "--dir", type=str, required=True, help="Directory containing chromosome-split MAF files.")
    parser.add_argument("-s", "--specie", type=str, required=True, help="Species name used in the MAF file.")

    args = parser.parse_args()

    # Validate input paths
    if not os.path.isfile(args.genePred):
        print(f"Error: genePred file not found: {args.genePred}")
        return
    if not os.path.isdir(args.dir):
        print(f"Error: Directory not found: {args.dir}")
        return

    # Parse genePred file and process entries
    try:
        gene_pred_data = parse_gene_pred(args.genePred)
        for entry in gene_pred_data:
            maf_parse(args.dir, entry, args.specie)
    except Exception as e:
        print(f"Error during processing: {e}")

if __name__ == "__main__":
    main()
