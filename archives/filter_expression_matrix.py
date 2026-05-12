import pandas as pd
import numpy as np
import argparse

def filter_by_expression(df, min_expr=1.0, min_samples=2):
    """
    过滤低表达基因：在至少 min_samples 个样本中表达值 >= min_expr。
    """
    expr_filter = (df >= min_expr).sum(axis=1) >= min_samples
    return df[expr_filter]

def filter_by_mad(df, mad_percentile=75):
    """
    手动计算每行 MAD，并保留 MAD 在前 mad_percentile% 的基因。
    """
    median = df.median(axis=1)
    mad = (df.sub(median, axis=0)).abs().median(axis=1)
    cutoff = np.percentile(mad, 100 - mad_percentile)
    return df[mad >= cutoff]


def main(args):
    # 读取表达矩阵
    df = pd.read_csv(args.input, sep='\t', index_col=0)
    
    # 过滤低表达
    df_filtered = filter_by_expression(df, min_expr=args.min_expr, min_samples=args.min_samples)
    
    # 保留 MAD 前百分比
    df_filtered = filter_by_mad(df_filtered, mad_percentile=args.mad_percentile)

    # 输出结果
    df_filtered.to_csv(args.output, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter low expression and low MAD genes from expression matrix")
    parser.add_argument("-i", "--input", required=True, help="输入基因表达矩阵文件 (TSV)")
    parser.add_argument("-o", "--output", required=True, help="输出文件名")
    parser.add_argument("--min_expr", type=float, default=1.0, help="表达值最低阈值")
    parser.add_argument("--min_samples", type=int, default=2, help="表达值超过阈值的最少样本数")
    parser.add_argument("--mad_percentile", type=float, default=75, help="保留 MAD 排名前多少百分比的基因")

    args = parser.parse_args()
    main(args)
