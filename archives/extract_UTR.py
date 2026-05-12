#!/usr/bin/
#########################################################################
# File Name: /home/chengyu/myscripts/extract_UTR.py
# Author: ChengYu
# Description: 
# Created Time: Tue 15 Oct 2024 04:28:27 PM CST
#########################################################################


import csv
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_gtf(gtf_file):
    transcripts = defaultdict(lambda: {'exons': [], 'cds': []})

    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue  # Skip comments
            fields = line.strip().split('\t')
            chrom, source, feature_type, start, end, score, strand, frame, attributes = fields
            start, end = int(start), int(end)
            
            # Extract transcript_id from attributes
            attr_dict = {key.strip(): value.strip('"') for key, value in 
                         (attr.split() for attr in attributes.split(';') if attr.strip())}
            transcript_id = attr_dict.get('transcript_id')

            if feature_type == 'exon':
                transcripts[transcript_id]['exons'].append((chrom, start, end, strand))
            elif feature_type == 'CDS':
                transcripts[transcript_id]['cds'].append((chrom, start, end, strand))

    return transcripts

def find_utr_regions(transcripts):
    utr_info = defaultdict(lambda: {'5UTR': [], '3UTR': []})

    for transcript_id, data in transcripts.items():
        exons = sorted(data['exons'], key=lambda x: x[1])
        cds = sorted(data['cds'], key=lambda x: x[1])

        if not cds or not exons:
            continue  # Skip transcripts with no CDS or exons

        cds_start = cds[0][1]
        cds_end = cds[-1][2]

        for chrom, exon_start, exon_end, strand in exons:
            if strand == "+":
                if exon_start < cds_start:  # Potential 5' UTR
                    utr_info[transcript_id]['5UTR'].append((chrom, exon_start, min(exon_end, cds_start - 1)))
                if exon_end > cds_end:  # Potential 3' UTR
                    utr_info[transcript_id]['3UTR'].append((chrom, max(exon_start, cds_end + 1), exon_end))
            elif strand == "-":  # On negative strand, reverse the logic
                if exon_end > cds_end:  # Potential 5' UTR
                    utr_info[transcript_id]['5UTR'].append((chrom, max(exon_start, cds_end + 1), exon_end))
                if exon_start < cds_start:  # Potential 3' UTR
                    utr_info[transcript_id]['3UTR'].append((chrom, exon_start, min(exon_end, cds_start - 1)))

    return utr_info

def extract_utr_sequences(utr_info, genome_fasta, output_fasta, output_csv):
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))

    with open(output_fasta, 'w') as fasta_out, open(output_csv, 'w', newline='') as csv_out:
        csv_writer = csv.writer(csv_out)
        csv_writer.writerow(['transcript_id', '5UTR_length', '3UTR_length'])

        for transcript_id, utrs in utr_info.items():
            utr_lengths = {'5UTR': 0, '3UTR': 0}

            for utr_type, regions in utrs.items():
                for chrom, start, end in regions:
                    seq = genome[chrom].seq[start-1:end]  # 1-based to 0-based
                    utr_lengths[utr_type] += len(seq)

                    fasta_out.write(f">{transcript_id}_{utr_type}\n{seq}\n")

            csv_writer.writerow([transcript_id, utr_lengths['5UTR'], utr_lengths['3UTR']])

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="Extract UTRs to fasta and length info to a csv file.")
    parser.add_argument('--gtf', required=True, help="GTF filepath")
    parser.add_argument('--genome', required=True, help="Genome FASTA filepath")
    parser.add_argument('--output_fasta', required=True, help="Output UTR FASTA filename.")
    parser.add_argument('--output_csv', required=True, help="Output UTR length CSV filename.")

    args = parser.parse_args()

    # 解析 GTF 文件
    transcripts = parse_gtf(args.gtf)
    
    # 查找 5' UTR 和 3' UTR 区域
    utr_info = find_utr_regions(transcripts)
    
    # 提取 UTR 序列并输出到文件
    extract_utr_sequences(utr_info, args.genome, args.output_fasta, args.output_csv)

if __name__ == "__main__":
    main()

