import os
import argparse
import pandas as pd
from collections import Counter

# 最大ASA参考表
MAX_ASA = {
    'A': 106.0, 'R': 248.0, 'N': 157.0, 'D': 163.0, 'C': 135.0,
    'Q': 198.0, 'E': 194.0, 'G': 84.0, 'H': 184.0, 'I': 169.0,
    'L': 164.0, 'K': 205.0, 'M': 188.0, 'F': 197.0, 'P': 136.0,
    'S': 130.0, 'T': 142.0, 'W': 227.0, 'Y': 222.0, 'V': 142.0
}

STRUCTURE_MAP = {
    'H': 'Helix', 'G': 'Helix', 'I': 'Helix',
    'E': 'Sheet', 'B': 'Sheet',
    'T': 'Turn', 'S': 'Turn'
}

def parse_dssp_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
    
    # 找到数据开始位置
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('  #  RESIDUE AA STRUCTURE'):
            data_start = i + 1
            break
    
    if data_start is None:
        print(f"[Warning] 无法找到数据头：{filepath}")
        return None

    structure_counts = {'Helix': 0, 'Sheet': 0, 'Turn': 0, 'Coil': 0}
    raw_structure_counter = Counter()
    relative_asa_list = []
    total_residues = 0

    for line in lines[data_start:]:
        if len(line) < 38:
            continue
        
        aa = line[13].strip()
        struct = line[16].strip()
        acc = line[34:38].strip()
        
        if not aa or aa == '!':
            continue
        
        max_asa = MAX_ASA.get(aa.upper())
        if not max_asa:
            continue
        
        try:
            acc = float(acc)
            relative_asa = acc / max_asa
            relative_asa_list.append(relative_asa)
        except ValueError:
            continue
        
        # 统计原始结构类型
        raw_structure_counter[struct] += 1
        
        # 统计归类后的结构类型
        struct_class = STRUCTURE_MAP.get(struct, 'Coil')
        structure_counts[struct_class] += 1
        total_residues += 1

    # 计算归类后结构的比例
    structure_ratios = {k: v / total_residues if total_residues else 0 for k, v in structure_counts.items()}
    
    # 计算平均相对ASA
    avg_relative_asa = sum(relative_asa_list) / len(relative_asa_list) if relative_asa_list else 0

    # 构建结果字典
    result = {
        'Filename': os.path.basename(filepath),
        'TotalResidues': total_residues,
        'AvgRelativeASA': avg_relative_asa,
    }
    
    # 添加归类后的结构统计（数量和比例）
    for k in ['Helix', 'Sheet', 'Turn', 'Coil']:
        result[f'{k}_Count'] = structure_counts[k]
        result[f'{k}_Ratio'] = structure_ratios[k]
    
    # 定义所有可能的原始结构类型
    all_raw_types = ['H', 'G', 'I', 'E', 'B', 'T', 'S', ' ', '-']
    
    # 添加原始结构统计（仅数量，不计算比例）
    for raw_type in all_raw_types:
        count = raw_structure_counter.get(raw_type, 0)
        result[f'Raw_{raw_type}_Count'] = count

    return result

def batch_process(folder_path, output_file):
    print(f"\n🔍 正在分析目录: {folder_path}\n")
    
    results_summary = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith('.dssp'):
            filepath = os.path.join(folder_path, filename)
            print(f"处理文件: {filename}")
            result = parse_dssp_file(filepath)
            if result:
                results_summary.append(result)
    
    if not results_summary:
        print("❌ 没有成功分析的文件。")
        return
    
    # 创建DataFrame并确保所有缺失值填充为0
    df = pd.DataFrame(results_summary)
    df = df.fillna(0)
    
    # 保存到CSV文件
    df.to_csv(output_file, index=False)
    
    print(f"\n✅ 所有统计完成：")
    print(f" - 合并结构摘要和原始结构输出: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='批量分析DSSP文件，统计结构类型、比例与ASA。',
        epilog='示例: python dssp_summary.py -i ./folder -o summary.csv'
    )
    
    parser.add_argument('-i', '--input', required=True, help='输入文件夹路径（含多个 .dssp 文件）')
    parser.add_argument('-o', '--output', default='dssp_summary.csv', help='输出文件路径（默认: dssp_summary.csv）')
    
    args = parser.parse_args()
    
    batch_process(args.input, args.output)

if __name__ == '__main__':
    main()
