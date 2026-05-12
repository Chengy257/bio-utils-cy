import argparse
import json

def parse_fasta(filename):
    """解析FASTA文件，返回包含完整header和序列的字典列表"""
    entries = []
    current_header = None
    current_sequence = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    entries.append({
                        'header': current_header,
                        'sequence': ''.join(current_sequence)
                    })
                    current_sequence = []
                current_header = line[1:].strip()  # 保留完整header（去除>符号和首尾空格）
            else:
                current_sequence.append(line)
        if current_header is not None:  # 处理最后一个条目
            entries.append({
                'header': current_header,
                'sequence': ''.join(current_sequence)
            })
    return entries

def main():
    parser = argparse.ArgumentParser(description='将FASTA文件转换为原始header的JSON格式')
    parser.add_argument('input_fasta', help='输入FASTA文件路径')
    parser.add_argument('output_json', help='输出JSON文件路径')
    args = parser.parse_args()

    # 解析FASTA文件
    fasta_entries = parse_fasta(args.input_fasta)

    # 构建JSON结构
    json_output = []
    for entry in fasta_entries:
        json_entry = {
            "name": entry['header'],  # 直接使用原始header
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": entry['sequence'],
                        "glycans": [],
                        "modifications": [],
                        "count": 1
                    }
                }
            ]
        }
        json_output.append(json_entry)

    # 写入JSON文件
    with open(args.output_json, 'w') as f:
        json.dump(json_output, f, indent=4, ensure_ascii=False)  # 保留非ASCII字符

    print(f"转换完成！已生成 {len(json_output)} 个条目到 {args.output_json}")

if __name__ == "__main__":
    main()
