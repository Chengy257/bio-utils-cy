import argparse
from concurrent.futures import ThreadPoolExecutor
import os

# 读取fasta文件，返回id->sequence字典
def read_fasta(fasta_file):
    seqs = {}
    with open(fasta_file) as f:
        current_id = None
        seq_list = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(seq_list).upper()
                current_id = line[1:].split()[0]
                seq_list = []
            else:
                seq_list.append(line)
        if current_id:
            seqs[current_id] = ''.join(seq_list).upper()
    return seqs

# 互补链转换
def rev_comp(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

# 解析genePred文件，返回每条转录本的基本信息
def parse_genepred_line(line):
    parts = line.strip().split('\t')
    chrom = parts[1]
    strand = parts[2]
    tx_start = int(parts[3])
    tx_end = int(parts[4])
    cds_start = int(parts[5])
    cds_end = int(parts[6])
    exon_count = int(parts[7])
    exon_starts = list(map(int, parts[8].strip(',').split(',')))
    exon_ends = list(map(int, parts[9].strip(',').split(',')))
    return {
        'line': line.strip(),
        'chrom': chrom,
        'strand': strand,
        'tx_start': tx_start,
        'tx_end': tx_end,
        'cds_start': cds_start,
        'cds_end': cds_end,
        'exon_count': exon_count,
        'exon_starts': exon_starts,
        'exon_ends': exon_ends
    }

# 根据exons从基因组序列中提取转录本序列
def get_transcript_seq(chrom_seq, record, strand):
    seq_fragments = []
    val_positions = []
    for start, end in zip(record['exon_starts'], record['exon_ends']):
        seq_fragments.append(chrom_seq[start:end])
        val_positions.extend(range(start, end))
    seq = ''.join(seq_fragments).upper()
    if strand == '-':
        seq = rev_comp(seq)
        val_positions = val_positions[::-1]
    return seq, val_positions

# 多线程处理函数：挖掘候选ORF并返回结果字符串，线程安全负责返回而非写文件
def process_transcript(record, chrom_seqs, start_codons, min_len):
    chrom = record['chrom']
    strand = record['strand']
    cds_start = record['cds_start']
    cds_end = record['cds_end']
    exon_starts = record['exon_starts']
    exon_ends = record['exon_ends']
    line_info = record['line']

    results_fa = []
    results_genepred = []

    if chrom not in chrom_seqs:
        # 缺失染色体序列，跳过
        return results_fa, results_genepred

    chrom_seq = chrom_seqs[chrom]

    transcript_seq, val_positions = get_transcript_seq(chrom_seq, record, strand)
    tlen = len(transcript_seq)

    # 由于Perl代码对于loc1和loc2定位基因CDS定义，稍加对应：
    # val_positions是基因组坐标数组，找cds_start和cds_end在vals中位置：
    loc1 = 0
    loc2 = 0
    for idx, pos in enumerate(val_positions):
        if pos == cds_start:
            loc1 = idx + 1  # 1-based
        if pos == cds_end - 1:
            loc2 = idx + 2  # perl代码这里是加2偏移，不改为保持一致

    # cds_start==cds_end时，loc1=loc2=0
    if cds_start == cds_end:
        loc1 = 0
        loc2 = 0

    ra = 0

    for i in range(tlen - 2):
        codon = transcript_seq[i:i+3]
        if codon in start_codons:
            j = i + 3
            while j <= tlen - 3:
                stop_codon = transcript_seq[j:j+3]
                if stop_codon in {'TAG', 'TAA', 'TGA'}:
                    orf_len = j - i + 3
                    if orf_len >= min_len:
                        seq_orf = transcript_seq[i:j+3]
                        ra += 1
                        loc3 = i + 1
                        loc4 = j + 4

                        # 判定ORF类型
                        if loc1 == 0:
                            orf_type = "noncoding"
                        elif loc3 == loc1 and loc4 == loc2:
                            orf_type = "canonical"
                        elif loc3 < loc1 and loc4 < loc1:
                            orf_type = "uORF"
                        elif loc3 < loc1 and loc4 >= loc1 and loc4 < loc2:
                            orf_type = "ouORF"
                        elif loc3 >= loc1 and loc4 < loc2:
                            orf_type = "iORF"
                        elif loc3 > loc1 and loc3 < loc2 and loc4 > loc2:
                            orf_type = "odORF"
                        elif loc3 > loc1 and loc4 == loc2:
                            orf_type = "truncation"
                        elif loc3 < loc1 and loc4 == loc2:
                            orf_type = "extension"
                        elif loc3 >= loc2:
                            orf_type = "dORF"
                        elif loc3 == loc1 and loc4 != loc2:
                            orf_type = "seqerror"
                        elif loc3 <= loc1 and loc4 > loc2:
                            orf_type = "readthrough"
                        else:
                            orf_type = "other"

                        # 生成id
                        id_str = f"{line_info.split()[0]}:{chrom}:{strand}|{ra}|{tlen}:{loc3}:{loc4}|{orf_type}|{codon}"

                        # fasta格式序列
                        fasta_entry = f">{id_str}\n{seq_orf}\n"
                        results_fa.append(fasta_entry)

                        # genePred格式输出，先拆分line_info
                        fields = line_info.split('\t')
                        if strand == '+':
                            start_out = val_positions[i]
                            end_out = val_positions[j+2] + 1
                        else:
                            start_out = val_positions[j+2]
                            end_out = val_positions[i] + 1

                        # 重新组合genePred字段，注意排列
                        # 之前部分字段是 parts 3~9，对应fields[3] ~ fields[9]
                        genepred_entry = (
                            f"{id_str}\t{chrom}\t{strand}\t"
                            f"{fields[3]}\t{fields[4]}\t"
                            f"{start_out}\t{end_out}\t"
                            f"{fields[7]}\t{fields[8]}\t{fields[9]}"
                        )
                        results_genepred.append(genepred_entry)

                    break  # 找到第一个终止密码子停止内层循环
                j += 3

    return results_fa, results_genepred

def main():
    parser = argparse.ArgumentParser(description="Annotate potential ORFs from genome and genePred transcriptome annotations.")
    parser.add_argument('-g', '--genome', required=True, help="Reference genome assembly file in fasta format.")
    parser.add_argument('-t', '--transcriptome', required=True, help="Reference transcriptome annotation file in genePred format.")
    parser.add_argument('-o', '--outdir', required=True, help="Output directory.")
    parser.add_argument('-s', '--startcodon', default="ATG/CTG/GTG/TTG/ACG", help="Start codons separated by '/', default: ATG/CTG/GTG/TTG/ACG")
    parser.add_argument('-l', '--minlen', type=int, default=6, help="Minimum candidate ORF length, default 6.")
    parser.add_argument('-p', '--threads', type=int, default=4, help="Number of threads, default 4.")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("Loading genome sequences...")
    chrom_seqs = read_fasta(args.genome)
    print(f"Loaded {len(chrom_seqs)} chromosomes/scaffolds.")

    start_codons = set(args.startcodon.split('/'))

    print("Loading transcriptome annotations...")
    transcripts = []
    with open(args.transcriptome) as f:
        for line in f:
            if line.strip():
                rec = parse_genepred_line(line)
                transcripts.append(rec)
    print(f"Loaded {len(transcripts)} transcripts.")

    print(f"Processing candidate ORFs with {args.threads} threads...")
    results_fa = []
    results_genepred = []

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_transcript, rec, chrom_seqs, start_codons, args.minlen) for rec in transcripts]
        for future in futures:
            fa_list, gp_list = future.result()
            results_fa.extend(fa_list)
            results_genepred.extend(gp_list)

    fa_out = os.path.join(args.outdir, 'candidateORF.fa')
    gp_out = os.path.join(args.outdir, 'candidateORF.genepred.txt')

    print(f"Writing {len(results_fa)} ORFs to {fa_out}...")
    with open(fa_out, 'w') as ffa:
        ffa.writelines(results_fa)

    print(f"Writing {len(results_genepred)} ORFs to {gp_out}...")
    with open(gp_out, 'w') as fgp:
        fgp.writelines(line + '\n' for line in results_genepred)

    print("Done!")

if __name__ == '__main__':
    main()
