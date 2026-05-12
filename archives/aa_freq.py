#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/aa_freq.py
# Author: ChengYu
# Description: 
# Created Time: Wed 12 Feb 2025 04:42:24 PM CST
#########################################################################
import argparse
import sys
from collections import defaultdict

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Calculate amino acid frequencies from a multi-FASTA file.',
        formatter_class=argparse.RawTextHelpFormatter
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
        '-c', '--ignore-case', action='store_true',
        help='Convert all amino acid characters to uppercase'
    )
    parser.add_argument(
        '-s', '--sort', choices=['alpha', 'freq'], default='alpha',
        help='Sort output by: alpha - alphabetical order (default)\n'
             '               freq - descending frequency order'
    )
    return parser.parse_args()

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

def main():
    args = parse_args()
    
    # 统计氨基酸频率
    counts = defaultdict(int)
    total = 0
    
    try:
        for seq in parse_fasta(args.input):
            processed_seq = seq.upper() if args.ignore_case else seq
            for aa in processed_seq:
                counts[aa] += 1
            total += len(processed_seq)
    except UnicodeDecodeError:
        sys.exit(f"Error: '{args.input}' is not a valid text file")
    
    if total == 0:
        sys.exit("Error: No valid sequences found in input file")
    
    # 计算频率
    frequencies = {aa: (count/total)*100 for aa, count in counts.items()}
    
    # 排序处理
    if args.sort == 'alpha':
        sorted_items = sorted(frequencies.items())
    else:
        sorted_items = sorted(frequencies.items(), 
                             key=lambda x: (-x[1], x[0]))
    
    # 写入输出文件
    try:
        with open(args.output, 'w') as f:
            f.write("Amino_acid\tCount\tFrequency(%)\n")
            for aa, freq in sorted_items:
                f.write(f"{aa}\t{counts[aa]}\t{freq:.2f}\n")
    except IOError:
        sys.exit(f"Error: Could not write to output file '{args.output}'")

    print(f"Successfully processed {total} amino acids")
    print(f"Results saved to: {args.output}")

if __name__ == '__main__':
    main()
