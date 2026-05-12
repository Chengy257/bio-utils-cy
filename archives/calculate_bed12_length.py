#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/calculate_bed12_length.py
# Author: ChengYu
# Description: 
# Created Time: Wed 08 Jan 2025 11:13:14 AM CST
#########################################################################

import argparse
from pybedtools import BedTool

def calculate_non_redundant_length(input_file):
    # 读取BED文件
    bed = BedTool(input_file)

    # 合并重叠区域
    merged_bed = bed.merge()

    # 计算总长度
    total_length = sum([int(end) - int(start) for start, end, _ in merged_bed])  # 只取前两列 start 和 end
    return total_length

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="计算BED文件中的非冗余基因组总长度")
    
    # 输入BED文件路径参数
    parser.add_argument('input_file', type=str, help="输入BED文件路径")
    
    # 解析命令行参数
    args = parser.parse_args()

    # 计算非冗余基因组总长度
    total_length = calculate_non_redundant_length(args.input_file)
    
    # 输出结果
    print(f"非冗余基因组总长度: {total_length}")

if __name__ == '__main__':
    main()



