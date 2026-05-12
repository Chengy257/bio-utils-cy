#!/bin/bash
#########################################################################
# File Name: batch_runDSSP.sh
# Author: ChengYu
# Description: 
# Created Time: Wed 11 Jun 2025 07:38:42 PM CST
#########################################################################

# 使用方式提示
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <PDB_DIR> <OUT_DIR> [NUM_CPU]"
    exit 1
fi

# 参数设置
PDB_DIR="$1"
OUT_DIR="$2"
NUM_CPU="${3:-20}"  # 默认使用 10 个线程
DSSP_BIN="/home/chengyu/soft/DSSP/dssp-4.5.3/build/mkdssp"
ENV_PATH="/home/chengyu/soft/miniconda3/envs/colabfold"

# 检查输入路径
if [ ! -d "$PDB_DIR" ]; then
    echo "Error: Input directory '$PDB_DIR' does not exist."
    exit 1
fi

# 创建输出路径
mkdir -p "$OUT_DIR"

# 激活 Conda 环境
source activate "$ENV_PATH"

# 检查 DSSP 可执行文件是否存在
if [ ! -x "$DSSP_BIN" ]; then
    echo "Error: DSSP binary not found or not executable at '$DSSP_BIN'"
    exit 1
fi

# 清空旧任务文件
BATCH_FILE="batch.jobs"
> "$BATCH_FILE"

# 遍历 PDB 文件
for PDB in "$PDB_DIR"/*.pdb; do
    [ -e "$PDB" ] || continue

    BASE_NAME="$(basename "$PDB")"
    OUT_FILE="$OUT_DIR/$BASE_NAME.dssp"

    # 如果输出文件已存在，则跳过
    if [ -s "$OUT_FILE" ]; then
        echo "Skipping existing: $OUT_FILE"
        continue
    fi

    # 修改 PDB 文件头（只插入一次）
    sed -i '1i\HEADER    ALPHAFOLD PREDICTION                      01-JAN-25   PRED\nCRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1' "$PDB"

    # 添加命令到任务列表
    printf "%s --calculate-accessibility --output-format dssp %s %s\n" \
        "$DSSP_BIN" "$PDB" "$OUT_FILE" >> "$BATCH_FILE"
done

# 并行执行任务
/opt/anaconda3/bin/ParaFly -c "$BATCH_FILE" -CPU "$NUM_CPU"



