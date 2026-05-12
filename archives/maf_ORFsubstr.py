#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author:       ChengYu
# Descriptions: Correctly extract CDS regions from a MAF file based on genePred data.
# Created:      2022/12/07
# Last Edited:  2024/11/21

import os
import argparse
from Bio.AlignIO import MafIO
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count



def calculate_cds_regions(cds_start, cds_end, exon_starts, exon_ends):
    """
    Correct CDS regions based on exon information.

    Args:
        cds_start (int): CDS start position from genePred.
        cds_end (int): CDS end position from genePred.
        exon_starts (list[int]): Start positions of exons.
        exon_ends (list[int]): End positions of exons.

    Returns:
        list[tuple]: List of corrected CDS regions as (start, end) tuples.
    """
    corrected_cds_regions = []
    for exon_start, exon_end in zip(exon_starts, exon_ends):
        # CDS region must intersect with exon region
        if exon_end > cds_start and exon_start < cds_end:
            # Calculate overlapping region
            corrected_start = max(cds_start, exon_start)
            corrected_end = min(cds_end, exon_end)
            corrected_cds_regions.append((corrected_start, corrected_end))
    return corrected_cds_regions


def maf_parse_task(args):
    """
    Task to extract and process CDS regions from chromosome-specific MAF files.

    Args:
        maf_dir (str): Directory containing chromosome-split MAF files.
        gene_pred (list): Parsed gene annotation data for a single entry.
    """
    maf_dir, gene_pred = args
    try:
        # Parse genePred data
        name, chr_name, strand, cds_start, cds_end, exon_count, exon_starts, exon_ends = (
            gene_pred[0], gene_pred[1], gene_pred[2],
            int(gene_pred[5]), int(gene_pred[6]),
            int(gene_pred[7]),
            list(map(int, filter(None, gene_pred[8].strip(",").split(",")))),
            list(map(int, filter(None, gene_pred[9].strip(",").split(","))))
        )
        strand = 1 if strand == "+" else -1

        # Calculate corrected CDS regions
        corrected_cds_regions = calculate_cds_regions(cds_start, cds_end, exon_starts, exon_ends)
        if not corrected_cds_regions:
            print(f"Warning: No valid CDS region found for {name}, skipping.")
            return

        # Prepare file paths
        maf_file = os.path.join(maf_dir, f"{chr_name}.maf")
        maf_index = f"{maf_file}.index"
        output_dna = f"{name}.dna.fa"
        output_protein = f"{name}.pep.fa"
        ref = chr_name

        if not os.path.isfile(maf_file):
            print(f"Warning: MAF file not found for chromosome {chr_name}, skipping entry {name}")
            return
        if not os.path.isfile(maf_index):
            print(f"Warning: Index file not found for {maf_file}, skipping entry {name}")
            return

        # Extract spliced CDS regions
        idx = MafIO.MafIndex(maf_index, maf_file, ref)
        start_positions = [region[0] for region in corrected_cds_regions]
        end_positions = [region[1] for region in corrected_cds_regions]
        try:
            region = idx.get_spliced(start_positions, end_positions, strand=strand)
        except KeyError:
            print(f"Warning: Target for indexing ({ref}) not found in MAF file {maf_file}, skipping entry {name}")
            return

        # Write extracted DNA FASTA
        AlignIO.write(region, output_dna, "fasta")

        # Translate to protein and write FASTA
        with open(output_dna, "r") as dna_in, open(output_protein, "w") as protein_out:
            for seq_record in SeqIO.parse(dna_in, "fasta"):
                protein_out.write(f">{seq_record.id}\n")
                # Remove gaps and translate
                protein_out.write(f"{Seq(str(seq_record.seq).replace('-', '')).translate()}\n")
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
    valid_entries = []
    try:
        with open(file_path, "r") as f:
            for line_number, line in enumerate(f, start=1):
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    print(f"Skipping invalid line {line_number}: {line.strip()}")
                    continue
                try:
                    int(fields[3])  # start position
                    int(fields[4])  # end position
                    int(fields[7])  # exon count
                except ValueError:
                    print(f"Skipping line with invalid numeric data at line {line_number}: {line.strip()}")
                    continue
                valid_entries.append(fields)
    except Exception as e:
        raise IOError(f"Error reading genePred file {file_path}: {e}")
    
    return valid_entries


def main():
    parser = argparse.ArgumentParser(description="Extract CDS regions from a MAF file based on genePred data.")
    parser.add_argument("-m", "--maf", type=str, required=True, help="Path to the whole MAF file.")
    parser.add_argument("-g", "--genePred", type=str, required=True, help="Path to the gene annotation file (genePred format).")
    parser.add_argument("-d", "--dir", type=str, required=True, help="Output directory for chromosome-split MAF files.")

    args = parser.parse_args()

    if not os.path.isfile(args.maf):
        print(f"Error: MAF file not found: {args.maf}")
        return
    if not os.path.isfile(args.genePred):
        print(f"Error: genePred file not found: {args.genePred}")
        return
    if not os.path.isdir(args.dir):
        os.makedirs(args.dir, exist_ok=True)

    try:
        gene_pred_data = parse_gene_pred(args.genePred)
        pool_args = [(args.dir, entry) for entry in gene_pred_data]

        num_workers = min(cpu_count(), len(pool_args))
        with Pool(processes=num_workers) as pool:
            pool.map(maf_parse_task, pool_args)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
