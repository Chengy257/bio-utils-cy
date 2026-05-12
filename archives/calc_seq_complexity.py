from Bio import SeqIO
import seq_entropies as se

def calc_complexity_for_seq(seq):
    # 将序列（字符串）切分为字符列表
    symbols = list(seq)
    # # LZ76 复杂度
    # _, lz76 = se.LZ76(symbols)
    # # LZ77 复杂度
    # _, lz77 = se.ZL77(symbols)
    # 块熵（以2-mer为例，base=2 比特）
    h2 = se.block_entropy(symbols, window=3, method="MLE", base=2)
    # 条件熵（同样以2-mer）
    c2 = se.block_cond_entropy(symbols, window=3, method="MLE", base=2)
    return h2, c2

def process_fasta(in_fasta, out_tsv):
    with open(out_tsv, 'w') as fout:
        # 输出标题行
        fout.write("seq_id\tBlockEntropy\tCondEntropy\n")
        for rec in SeqIO.parse(in_fasta, "fasta"):
            h2, c2 = calc_complexity_for_seq(str(rec.seq))
            fout.write(f"{rec.id}\t{h2:.4f}\t{c2:.4f}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python calc_seq_complexity.py <input.fasta> <output.tsv>")
        sys.exit(1)
    process_fasta(sys.argv[1], sys.argv[2])
