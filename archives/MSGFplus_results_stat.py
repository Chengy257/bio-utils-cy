#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/MSGFplus_results_stat.py
# Author: ChengYu
# Description: 
# Created Time: Wed 19 Feb 2025 05:39:53 PM CST
#########################################################################
import argparse
import pandas as pd
import os

def process_protein_data(args):
    """
    处理蛋白质数据并生成统计报告
    """
    try:
        # 读取输入文件
        df = pd.read_csv(args.input, sep='\t')
        print(f"成功读取输入文件: {args.input} (共{len(df)}行)")
    except FileNotFoundError:
        raise SystemExit(f"错误: 输入文件 {args.input} 不存在")
    except Exception as e:
        raise SystemExit(f"读取文件失败: {str(e)}")

    # 列名校验
    required_columns = ['Protein', 'Peptide', 'SpecId']
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        raise SystemExit(f"缺少必要列: {', '.join(missing_cols)}")

    # FDR过滤
    if args.fdr is not None:
        if 'QValue' not in df.columns:
            raise SystemExit("未找到QValue列，无法执行FDR过滤")
        original_count = len(df)
        df = df[df['PepQValue'] < args.fdr]
        print(f"应用FDR过滤 (<{args.fdr})，保留{len(df)}行 (过滤掉{original_count - len(df)}行)")

    # 处理多蛋白归属
    if args.split_proteins:
        df = df.assign(Protein=df['Protein'].str.split(';')).explode('Protein')
        print("已拆分多归属蛋白质条目")

    # 核心统计逻辑
    try:
        # 预处理数据
        df = df[required_columns].drop_duplicates()

        # 计算肽段-蛋白映射数
        peptide_counts = df.groupby('Peptide')['Protein'].nunique().reset_index()
        peptide_counts.columns = ['Peptide', 'ProteinCount']

        # 合并数据
        merged_df = pd.merge(df, peptide_counts, on='Peptide')

        # 定义统计函数
        def get_protein_stats(group):
            return pd.Series({
                'SpectraCount': group['SpecId'].nunique(),
                'TotalPeptides': group['Peptide'].nunique(),
                'UniquePeptides': group[group['ProteinCount'] == 1]['Peptide'].nunique()
            })

        # 生成统计结果
        protein_stats = merged_df.groupby('Protein').apply(get_protein_stats).reset_index()
    except Exception as e:
        raise SystemExit(f"数据处理失败: {str(e)}")

    # 保存结果
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    protein_stats.to_csv(args.output, index=False)
    print(f"成功生成报告: {args.output}")
    print(f"统计蛋白数量: {len(protein_stats)}")

def main():
    # 配置命令行参数
    parser = argparse.ArgumentParser(
        description="MSGF+蛋白质组学数据统计工具",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help="输入文件路径 (TSV格式)")
    parser.add_argument('-o', '--output', default='protein_report.csv',
                        help="输出文件路径 (默认: protein_report.csv)")
    parser.add_argument('--fdr', type=float, 
                        help="FDR过滤阈值 (例如: 0.01)")
    parser.add_argument('--split-proteins', action='store_true',
                        help="拆分多归属蛋白条目 (用分号分隔的蛋白列表)")
    parser.add_argument('--show-example', action='store_true',
                        help="显示示例输入格式")

    # 示例显示处理
    args = parser.parse_args()
    if args.show_example:
        print("示例输入列：")
        print("Protein\tPeptide\tSpecId\tQValue")
        print("sp|P12345\tAASLLK\tspectrum_001\t0.005")
        print("sp|P12345;sp|Q67890\tLLKASD\tspectrum_002\t0.02")
        return

    # 执行处理流程
    process_protein_data(args)

if __name__ == "__main__":
    main()
