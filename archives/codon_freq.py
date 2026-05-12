#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/codon_freq.py
# Author: ChengYu
# Description: 
# Created Time: Wed 12 Feb 2025 10:28:55 PM CST
#########################################################################
import argparse
import sys
from collections import defaultdict

VALID_NUCLEOTIDES = {'A', 'T', 'C', 'G', 'U'}

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Calculate codon frequencies from a multi-FASTA nucleotide file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Input FASTA file path (required)'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output file path (required)'
    )
    parser.add_argument(
        '-s', '--sort', choices=['alpha', 'freq'], default='alpha',
        help='Sort output by: alpha - alphabetical order, freq - descending frequency'
    )
    parser.add_argument(
        '-t', '--truncate', action='store_true',
        help='Truncate sequences to length divisible by 3 instead of raising error'
    )
    parser.add_argument(
        '--skip-invalid', action='store_true',
        help='Skip codons containing invalid nucleotides instead of raising error'
    )
    parser.add_argument(
        '-u', '--uppercase', action='store_true',
        help='Convert all nucleotides to uppercase'
    )
    return parser.parse_args()

def validate_codon(codon):
    """验证密码子是否有效"""
    return all(n in VALID_NUCLEOTIDES for n in codon)

def parse_fasta(file_path):
    """解析FASTA文件，生成序列迭代器"""
    try:
        with open(file_path, 'r') as f:
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        yield ''.join(current_seq)
                        current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                yield ''.join(current_seq)
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{file_path}' not found")

def process_sequence(seq, args):
    """处理单个序列，返回有效密码子列表"""
    seq = seq.upper() if args.uppercase else seq
    seq_len = len(seq)
    
    # 处理长度不是3的倍数的情况
    if seq_len % 3 != 0:
        if args.truncate:
            seq = seq[:seq_len - (seq_len % 3)]
        else:
            sys.exit(f"Error: Sequence length {seq_len} is not divisible by 3")

    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    
    # 过滤无效密码子
    valid_codons = []
    for codon in codons:
        if validate_codon(codon):
            valid_codons.append(codon)
        elif not args.skip_invalid:
            sys.exit(f"Error: Invalid codon '{codon}' found")
    
    return valid_codons

def main():
    args = parse_args()
    codon_counts = defaultdict(int)
    total = 0

    try:
        for seq in parse_fasta(args.input):
            codons = process_sequence(seq, args)
            for codon in codons:
                codon_counts[codon] += 1
            total += len(codons)
    except UnicodeDecodeError:
        sys.exit(f"Error: '{args.input}' is not a valid text file")

    if total == 0:
        sys.exit("Error: No valid codons found in input file")

    # 计算频率
    frequencies = {codon: (count/total)*100 for codon, count in codon_counts.items()}

    # 排序处理
    if args.sort == 'alpha':
        sorted_items = sorted(frequencies.items())
    else:
        sorted_items = sorted(frequencies.items(), 
                            key=lambda x: (-x[1], x[0]))

    # 写入输出文件
    try:
        with open(args.output, 'w') as f:
            f.write("Codon\tCount\tFrequency(%)\n")
            for codon, freq in sorted_items:
                f.write(f"{codon}\t{codon_counts[codon]}\t{freq:.4f}\n")
    except IOError:
        sys.exit(f"Error: Could not write to output file '{args.output}'")

    print(f"Successfully processed {total} codons")
    print(f"Results saved to: {args.output}")

if __name__ == '__main__':
    main()
