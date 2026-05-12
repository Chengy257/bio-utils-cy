from Bio import SeqIO
from CAI import CAI, relative_adaptiveness
import sys

def filter_cds_sequences(fasta_path, filtered_path, seq_type="reference"):
    """
    过滤非3倍数长度的序列，写入filtered_path
    """
    count_total = 0
    count_filtered = 0
    with open(filtered_path, "w") as fout:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            count_total += 1
            seq_len = len(rec.seq)
            if seq_len % 3 == 0:
                fout.write(f">{rec.id}\n{str(rec.seq)}\n")
                count_filtered += 1
            else:
                print(f"Skipped {seq_type} seq {rec.id} length={seq_len} (not divisible by 3)")
    print(f"{seq_type.capitalize()} sequences: total={count_total}, kept={count_filtered}, skipped={count_total - count_filtered}")
    return filtered_path

def is_valid_cds(seq):
    """
    判断序列是否长度是3的倍数且只含有标准碱基ATGC
    """
    bases = set("ATGC")
    seq = seq.upper()
    if len(seq) % 3 != 0:
        return False
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if any(b not in bases for b in codon):
            return False
    return True

def filter_target_sequences(input_fasta, filtered_path):
    """
    过滤目标序列：长度和碱基均过滤，写入filtered_path
    """
    count_total = 0
    count_filtered = 0
    with open(filtered_path, "w") as fout:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            count_total += 1
            seq = str(rec.seq).upper()
            if is_valid_cds(seq):
                fout.write(f">{rec.id}\n{seq}\n")
                count_filtered += 1
            else:
                print(f"Skipped target seq {rec.id} due to invalid CDS (non-ATGC or length not multiple of 3)")
    print(f"Target sequences: total={count_total}, kept={count_filtered}, skipped={count_total - count_filtered}")
    return filtered_path

def batch_cai_with_filter(input_fasta, ref_fasta, output_tsv):
    # 过滤参考序列
    filtered_ref = "filtered_reference.fasta"
    filter_cds_sequences(ref_fasta, filtered_ref, seq_type="reference")

    # 过滤目标序列
    filtered_input = "filtered_input.fasta"
    filter_target_sequences(input_fasta, filtered_input)

    # 读取过滤后的参考序列，计算权重
    reference = [str(rec.seq).upper() for rec in SeqIO.parse(filtered_ref, "fasta")]
    weights = relative_adaptiveness(sequences=reference)

    # 计算并输出CAI
    with open(output_tsv, "w") as fout:
        fout.write("seq_id\tCAI\n")
        for rec in SeqIO.parse(filtered_input, "fasta"):
            seq = str(rec.seq).upper()
            try:
                cai_value = CAI(seq, weights=weights)
            except Exception as e:
                print(f"Error computing CAI for {rec.id}: {e}")
                continue
            fout.write(f"{rec.id}\t{cai_value:.6f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calc_cai_batch_robust.py <input.fasta> <reference.fasta> <output.tsv>")
        sys.exit(1)
    batch_cai_with_filter(sys.argv[1], sys.argv[2], sys.argv[3])
