#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/ExtractBedRegionsFromMAF.py
# Author: ChengYu
# Description: 
# Created Time: Sat 28 Dec 2024 04:32:19 PM CST
#########################################################################

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bx.align import maf

def read_bed12(bed_file):
    """解析bed12文件，返回包含CDS信息的字典"""
    cds_regions = {}
    with open(bed_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom, start, end, name, _, strand, _, _, _, block_count, block_sizes, block_starts = fields[:12]
            block_sizes = [int(size) for size in block_sizes.split(",") if size]
            block_starts = [int(start) for start in block_starts.split(",") if start]

            cds_regions[name] = {
                "chrom": chrom,
                "strand": strand,
                "blocks": [
                    (int(start) + start_offset, int(start) + start_offset + size)
                    for start_offset, size in zip(block_starts, block_sizes)
                ]
            }
    # print(cds_regions)
    return cds_regions

def extract_cds_from_maf(maf_file, cds_regions, reference_species):
    """从MAF文件中提取对应的CDS区域序列"""
    cds_sequences = {}

    with open(maf_file) as maf_f:
        maf_reader = maf.Reader(maf_f)

        for alignment in maf_reader:
            for name, info in cds_regions.items():
                ref_chrom = f"{reference_species}.{info['chrom']}"
                for comp in alignment.components:
                    # 匹配参考物种和染色体
                    if comp.src.startswith(ref_chrom):
                        print(f"Processing {comp.src} for {name}")  # 调试输出
                        cds_seq = ""
                        for block_start, block_end in info["blocks"]:
                            try:
                                # 获取序列起止位置
                                chrom, coords = comp.src.split(":")
                                start, size = map(int, coords.split("-"))
                                print(f"Block: {block_start}-{block_end}, Alignment: {start}-{start+size}")  # 调试输出

                                if not (block_end < start or block_start > start + size):
                                    aln_start = max(block_start, start)
                                    aln_end = min(block_end, start + size)
                                    seq_start = aln_start - start
                                    seq_end = aln_end - start
                                    cds_seq += comp.text[seq_start:seq_end]
                            except ValueError as e:
                                print(f"Error processing {comp.src}: {e}")
                                continue

                        if cds_seq:  # 如果有序列，处理负链
                            if info["strand"] == "-":
                                cds_seq = str(Seq(cds_seq).reverse_complement())

                            cds_sequences[name] = cds_seq.replace("-", "")  # 去掉gap
                        else:
                            print(f"No sequence extracted for {name}")
    return cds_sequences

def write_fasta(file_name, sequences, translate=False):
    """将序列写入fasta文件"""
    records = []
    for name, seq in sequences.items():
        if translate:
            seq = str(Seq(seq).translate(to_stop=True))
        records.append(SeqRecord(Seq(seq), id=name, description=""))
    SeqIO.write(records, file_name, "fasta")

def main():
    parser = argparse.ArgumentParser(description="从MAF文件中提取CDS并翻译")
    parser.add_argument("--bed", required=True, help="输入的BED12文件路径")
    parser.add_argument("--maf", required=True, help="输入的MAF文件路径")
    parser.add_argument("--reference", required=True, help="参考物种名，例如'hg38'")
    parser.add_argument("--cds_fasta", required=True, help="输出的CDS fasta文件路径")
    parser.add_argument("--protein_fasta", required=True, help="输出的蛋白质序列fasta文件路径")

    args = parser.parse_args()

    # 解析bed12文件
    cds_regions = read_bed12(args.bed)

    # 从MAF文件中提取CDS序列
    cds_sequences = extract_cds_from_maf(args.maf, cds_regions, args.reference)

    # 写入CDS fasta文件
    write_fasta(args.cds_fasta, cds_sequences)

    # 写入蛋白质序列fasta文件
    write_fasta(args.protein_fasta, cds_sequences, translate=True)

    print("CDS和蛋白质序列已生成。")

if __name__ == "__main__":
    main()
