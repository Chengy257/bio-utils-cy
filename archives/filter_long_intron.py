#!/usr/bin/env python3
"""
Filter transcripts with abnormally long introns from a GTF/GFF file.

Usage:
    python filter_long_intron.py -i input.gtf -o output.gtf -m 20000

Author: (optimized by ChatGPT)
"""

import argparse
import sys
import os
from tqdm import tqdm
from BCBio import GFF

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter transcripts with introns longer than a threshold."
    )
    parser.add_argument("-i", "--input", required=True, help="Input GTF/GFF file")
    parser.add_argument("-o", "--output", required=True, help="Output filtered GTF file")
    parser.add_argument("-m", "--max_intron", type=int, default=20000,
                        help="Maximum intron length (bp) [default: 20000]")
    return parser.parse_args()

def filter_gtf(input_file, output_file, max_intron):
    if not os.path.exists(input_file):
        sys.exit(f"❌ Error: Input file '{input_file}' not found.")

    print(f"🔍 Filtering {input_file} ... (max intron = {max_intron:,} bp)")

    in_handle = open(input_file)
    out_handle = open(output_file, "w")

    total_tx, filtered_tx = 0, 0

    # use tqdm for progress bar
    for rec in tqdm(GFF.parse(in_handle), desc="Processing sequences", unit="seq"):
        new_feats = []
        for feat in rec.features:
            total_tx += 1
            exons = [f.location for f in feat.sub_features if f.type == "exon"]
            if len(exons) > 1:
                intron_lengths = [exons[i+1].start - exons[i].end for i in range(len(exons)-1)]
                if any(l > max_intron for l in intron_lengths):
                    filtered_tx += 1
                    continue  # skip this transcript
            new_feats.append(feat)
        rec.features = new_feats
        GFF.write([rec], out_handle)

    in_handle.close()
    out_handle.close()

    print("\n📊 Filtering Summary")
    print("--------------------")
    print(f"Total transcripts:   {total_tx:,}")
    print(f"Filtered (too long): {filtered_tx:,}")
    print(f"Retained:            {total_tx - filtered_tx:,}")
    print(f"✅ Output written to: {output_file}")

def main():
    args = parse_args()
    filter_gtf(args.input, args.output, args.max_intron)

if __name__ == "__main__":
    main()
