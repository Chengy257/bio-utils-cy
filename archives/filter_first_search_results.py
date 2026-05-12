#!/usr/bin/
#########################################################################
# File Name: filter_first_search_results.py
# Author: ChengYu
# Description: 
# Created Time: Fri 07 Feb 2025 04:55:11 PM CST
#########################################################################

import argparse
import pandas as pd
import pyteomics.mgf

def extract_scan_from_specid(specid):
    """
    从 SpecId 中提取 scan 的数字部分（等号后的数字）。
    """
    try:
        return int(specid.split('=')[0])
    except IndexError:
        raise ValueError(f"Invalid SpecId format: {specid}")

def filter_first_search_results(results, mgf, filtered_mgf, threshold_score=0.01):
    """
    过滤第一次搜索结果，筛选出得分低于阈值的谱图并写入新的 MGF 文件。
    
    参数:
    results (str): 输入的第一次搜索结果 TSV 文件路径。
    mgf (str): 输入的原始 MGF 文件路径。
    filtered_mgf (str): 输出的过滤后的 MGF 文件路径。
    threshold_score (float): 过滤的阈值，默认是 0.01。
    """
    # 读取第一次搜库的结果
    df = pd.read_csv(results, sep="\t")
    
    # 筛选得分低于阈值的谱图
    filtered_df = df[df['PepQValue'] <= threshold_score]
    
    # 提取所有符合条件的谱图 SpecId，并提取出其中的 scan 值
    filtered_spectra = {extract_scan_from_specid(specid) for specid in filtered_df['ScanNum']}

    # 读取原始 MGF 文件并筛选谱图
    with open(mgf, 'r') as f, open(filtered_mgf, 'w') as out_f:
        for spectrum in pyteomics.mgf.read(f):
            # 提取 TITLE 行中的 scan 值
            title = spectrum.get('params', {}).get('title', '')
            scan_number = None
            if title.startswith('scan='):
                try:
                    scan_number = int(title.split('=')[1])
                except ValueError:
                    pass  # 如果解析出错，跳过该谱图
            # 如果 scan_number 与过滤结果中的值匹配，跳过该谱图
            if scan_number not in filtered_spectra:
                pyteomics.mgf.write(spectrum, out_f)

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Filter first search results and generate a filtered MGF file.")
    parser.add_argument('results', type=str, help="Path to the first search results TSV file.")
    parser.add_argument('mgf', type=str, help="Path to the original MGF file.")
    parser.add_argument('filtered_mgf', type=str, help="Path to the filtered MGF output file.")
    parser.add_argument('--threshold', type=float, default=0.01, help="Threshold score for filtering. Default is 0.01.")

    args = parser.parse_args()
    
    # 调用过滤函数
    filter_first_search_results(args.results, args.mgf, args.filtered_mgf, args.threshold)

if __name__ == '__main__':
    main()
