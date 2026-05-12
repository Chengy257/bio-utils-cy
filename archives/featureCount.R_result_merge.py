#!/usr/bin/env python3
#########################################################################
# File Name: /home/chengyu/myscripts/featureCount.R_result_merge.py
# Author: ChengYu
# Description: 
# Created Time: Fri 14 Nov 2025 11:01:52 AM CST
#########################################################################

import os
import sys
import pandas as pd
from pathlib import Path

def read_count_file(path):
    df = pd.read_csv(path, sep="\t")
    df = df[["id", "effLength", "counts", "fpkm", "tpm"]]
    return df

def merge_metric(count_files, metric, output_file):
    """Merge counts / fpkm / tpm into desired output file."""
    merged = None

    for f in count_files:
        sample = Path(f).stem
        df = read_count_file(f)[["id", metric]].rename(columns={metric: sample})

        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on="id", how="outer")

    merged.to_csv( output_file, sep="\t", index=False)
    print(f"[OK] {output_file}")

# def merge_log_files(log_files, output_file):
#     """Merge Assigned counts from each log file."""
#     rows = []

#     for f in log_files:
#         sample = Path(f).stem
#         assigned = None

#         with open(f) as fh:
#             for line in fh:
#                 if line.startswith("Assigned"):
#                     assigned = line.strip().split("\t")[1]
#                     break

#         if assigned is not None:
#             rows.append([sample, assigned])

#     df = pd.DataFrame(rows, columns=["sample", "Assigned"])
#     df.to_csv(output_file, sep="\t", index=False)
#     print(f"[OK] {output_file}")

def merge_log_files(log_files, output_file):
    """Merge all lines from *.log files into a single table."""
    dfs = []
    for f in log_files:
        sample = Path(f).stem
        df = pd.read_csv(f, sep="\t", header=None, names=["metric", sample])
        df.set_index("metric", inplace=True)
        dfs.append(df)

    # 按行合并所有样本
    merged = pd.concat(dfs, axis=1)
    merged.to_csv(output_file, sep="\t")
    print(f"[OK] {output_file} generated.")



def main(input_dir):
    input_dir = Path(input_dir)

    count_files = sorted(input_dir.glob("*.count"))
    log_files = sorted(input_dir.glob("*.log"))

    if not count_files:
        print("No *.count files found.")
        sys.exit(1)

    print(f"Found {len(count_files)} count files.")

    # ------------------------------
    # effLength：直接取第一个文件
    # ------------------------------
    df0 = read_count_file(count_files[0])[["id", "effLength"]]
    df0.to_csv(input_dir / "effLength.txt", sep="\t", index=False)
    print("[OK] effLength.txt generated.")

    # ------------------------------
    # 三个最终合并文件名（完全符合你要求）
    # ------------------------------
    merge_metric(count_files, "counts", input_dir / "count.matrix.tsv")
    merge_metric(count_files, "fpkm",   input_dir / "GeneExpression_FPKM.xls")
    merge_metric(count_files, "tpm",    input_dir / "GeneExpression_TPM.xls")

    # ------------------------------
    # log 文件（Assigned）
    # ------------------------------
    if log_files:
        merge_log_files(log_files, input_dir / "GeneCount_Assigned_logs.xls")

    print("All done!")

if __name__ == "__main__":
    main(sys.argv[1])
