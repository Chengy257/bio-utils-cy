#!/usr/bin/env python3
#########################################################################
# File Name: /home/chengyu/myscripts/get5UTR_genePred2BED.py
# Author: ChengYu
# Description: 
# Created Time: Wed 11 Dec 2024 05:46:44 PM CST
#########################################################################

import argparse
from multiprocessing import Pool

def parse_genePred_line(line):
    """Parse a single line of genePred and extract necessary information."""
    fields = line.strip().split("\t")
    gene_id = fields[0]
    chrom = fields[1]
    strand = fields[2]
    tx_start = int(fields[3])
    tx_end = int(fields[4])
    cds_start = int(fields[5])
    cds_end = int(fields[6])
    exon_starts = list(map(int, fields[8].strip(",").split(",")))
    exon_ends = list(map(int, fields[9].strip(",").split(",")))
    return gene_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends

def extract_utr_bed12(gene_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends, utr_type):
    """Extract UTR regions and format them as BED12."""
    utr_regions = []
    block_sizes = []
    block_starts = []
    chrom_start = None
    chrom_end = None

    if utr_type == "5UTR":
        for start, end in zip(exon_starts, exon_ends):
            if end <= cds_start:
                if chrom_start is None:
                    chrom_start = start
                chrom_end = end
                block_sizes.append(end - start)
                block_starts.append(start - chrom_start)
            elif start < cds_start:
                if chrom_start is None:
                    chrom_start = start
                chrom_end = cds_start
                block_sizes.append(cds_start - start)
                block_starts.append(start - chrom_start)
                break
            else:
                break
    elif utr_type == "3UTR":
        for start, end in zip(exon_starts, exon_ends):
            if start >= cds_end:
                if chrom_start is None:
                    chrom_start = start
                chrom_end = end
                block_sizes.append(end - start)
                block_starts.append(start - chrom_start)
            elif end > cds_end:
                if chrom_start is None:
                    chrom_start = cds_end
                chrom_end = end
                block_sizes.append(end - cds_end)
                block_starts.append(cds_end - chrom_start)

    if chrom_start is None or chrom_end is None:
        return None

    return [
        chrom, chrom_start, chrom_end, gene_id, 0, strand,
        chrom_start, chrom_start, 0,
        len(block_sizes), ",".join(map(str, block_sizes)), ",".join(map(str, block_starts))
    ]

def process_genePred_line(line, utr_type):
    """Process a single genePred line and extract the specified UTR."""
    gene_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends = parse_genePred_line(line)
    return extract_utr_bed12(gene_id, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_starts, exon_ends, utr_type)

def genePred_to_bed12(input_file, output_file, utr_type, threads):
    """Convert genePred file to BED12 file with UTR regions."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()

        # Process lines in parallel if threads > 1
        if threads > 1:
            with Pool(threads) as pool:
                results = pool.starmap(process_genePred_line, [(line, utr_type) for line in lines])
        else:
            results = [process_genePred_line(line, utr_type) for line in lines]

        for result in results:
            if result:
                outfile.write("\t".join(map(str, result)) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Convert genePred file to BED12 file with UTR regions.")
    parser.add_argument("-i", "--input", required=True, help="Input genePred file")
    parser.add_argument("-o", "--output", required=True, help="Output BED12 file")
    parser.add_argument("--utr", choices=["5UTR", "3UTR", "both"], default="5UTR", 
                        help="Type of UTR to extract (default: 5UTR)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use (default: 1)")

    args = parser.parse_args()

    try:
        if args.utr == "both":
            # Handle both 5'UTR and 3'UTR cases separately
            output_5utr = args.output.replace(".bed12", "_5UTR.bed12")
            output_3utr = args.output.replace(".bed12", "_3UTR.bed12")
            print(f"Extracting 5'UTR to {output_5utr}")
            genePred_to_bed12(args.input, output_5utr, "5UTR", args.threads)
            print(f"Extracting 3'UTR to {output_3utr}")
            genePred_to_bed12(args.input, output_3utr, "3UTR", args.threads)
        else:
            # Handle a single UTR type
            print(f"Extracting {args.utr} to {args.output}")
            genePred_to_bed12(args.input, args.output, args.utr, args.threads)
        print("Conversion completed.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
