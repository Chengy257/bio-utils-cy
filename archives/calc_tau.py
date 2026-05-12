#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

def calc_tau(values):
    values = np.array(values, dtype=float)
    if np.all(values == 0):
        return np.nan  # 全0表达返回NaN
    max_expr = np.max(values)
    tau = np.sum(1 - (values / max_expr)) / (len(values) - 1)
    return tau

def main(input_file, output_file):
    # 读取输入文件（制表符分隔）
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.rename(columns={0: 'Gene'}, inplace=True)
    
    # 计算 tau 值
    df['tau'] = df.iloc[:, 1:].apply(calc_tau, axis=1)
    
    # 输出结果
    df_out = df[['Gene', 'tau']]
    df_out.to_csv(output_file, sep='\t', index=False)
    print(f"✅ 结果已保存到: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python calc_tau.py <expression_matrix.txt> <output_tau.txt>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

