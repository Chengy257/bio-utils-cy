#!/usr/bin/python
#########################################################################
# File Name: maf_ORFsubstr.2412.py
# Author: ChengYu
# Description: 
# Created Time: Sat 28 Dec 2024 10:26:33 PM CST
#########################################################################
import os
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


# 初始化日志
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

def parse_bed12(bed12_file):
    """
    Parse a BED12 file into a list of entries.

    Args:
        bed12_file (str): Path to the BED12 file.

    Returns:
        list: Parsed BED12 data entries.
    """
    try:
        with open(bed12_file, "r") as f:
            entries = []
            for line in f:
                if not line.strip():
                    continue
                fields = line.strip().split("\t")
                chr_name, start, end, name, strand, block_sizes, block_starts = (
                    fields[0], int(fields[1]), int(fields[2]), fields[3],
                    fields[5], fields[10], fields[11]
                )
                block_sizes = list(map(int, block_sizes.strip(",").split(",")))
                block_starts = list(map(int, block_starts.strip(",").split(",")))
                entries.append({
                    "chr": chr_name,
                    "start": start,
                    "end": end,
                    "name": name,
                    "strand": 1 if strand == "+" else -1,
                    "block_sizes": block_sizes,
                    "block_starts": block_starts,
                })
            return entries
    except Exception as e:
        raise IOError(f"Error reading BED12 file {bed12_file}: {e}")

def filter_species(region, ref):
    """
    Filter out species with no aligned bases in the MAF region.

    Args:
        region (MultipleSeqAlignment): MAF alignment region.
        ref (str): Reference species identifier.

    Returns:
        MultipleSeqAlignment: Filtered alignment region with only relevant species.
    """
    from Bio.Align import MultipleSeqAlignment
    filtered_records = []
    for record in region:
        if str(record.seq).replace("-", "").strip():  # 检查是否存在比对碱基
            filtered_records.append(record)
    return MultipleSeqAlignment(filtered_records)

def maf_parse(maf_dir, entry, species):
    """
    Extract and process block regions from a MAF file based on BED12 data,
    and output the results to both FASTA and MAF format.

    Args:
        maf_dir (str): Directory containing chromosome-split MAF files.
        entry (dict): BED12 data for a single entry.
        species (str): Species name used in the MAF file.
    """
    try:
        # 从BED12条目中提取基本信息
        name = entry["name"]
        chr_name = entry["chr"]
        strand = entry["strand"]
        block_sizes = entry["block_sizes"]
        block_starts = entry["block_starts"]
        maf_file = os.path.join(maf_dir, f"{chr_name}.maf")
        maf_index = f"{maf_file}.index"
        output_dna = f"{name}.dna.fa"
        output_protein = f"{name}.pep.fa"
        output_maf = f"{name}.maf"
        ref = f"{species}.{chr_name}"

        if not os.path.isfile(maf_file):
            logging.warning(f"MAF file not found: {maf_file}")
            return

        # 创建MafIndex对象
        idx = MafIO.MafIndex(maf_index, maf_file, ref)

        # 计算每个block的绝对区间
        start_positions = [entry["start"] + start for start in block_starts]
        end_positions = [start + size for start, size in zip(start_positions, block_sizes)]

        try:
            # 使用 get_spliced 方法提取拼接区间
            region = idx.get_spliced(start_positions, end_positions, strand=strand)
        except Exception as e:
            logging.error(f"Error extracting spliced region for {name}: {e}")
            return

        # 检查是否成功提取区域
        if not region:
            logging.warning(f"No valid blocks extracted for {name}. Skipping.")
            return

        # 将提取的拼接区域写入到 MAF 文件
        with open(output_maf, "w") as maf_out:
            AlignIO.write(region, maf_out, "maf")
        logging.info(f"MAF output written for {name}: {output_maf}")

        # 移除无比对的物种
        filtered_region = filter_species(region, ref)
        if not filtered_region:
            logging.warning(f"No relevant sequences for {name}. Skipping.")
            return

        # 写入DNA序列到FASTA
        AlignIO.write(filtered_region, output_dna, "fasta")

        # 翻译DNA为蛋白质
        with open(output_dna, "r") as dna_in, open(output_protein, "w") as protein_out:
            for seq_record in SeqIO.parse(dna_in, "fasta"):
                seq_ungapped = seq_record.seq.replace("-", "")
                protein_out.write(f">{seq_record.id}\n{seq_ungapped.translate()}\n")
        logging.info(f"Processed BED entry: {name}")
    except Exception as e:
        logging.error(f"Error processing BED entry {entry['name']}: {e}")


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract subregions from MAF files and output as FASTA.")
    parser.add_argument("-b", "--bed12", type=str, required=True, help="Path to the BED12 file.")
    parser.add_argument("-d", "--dir", type=str, required=True, help="Directory containing chromosome-split MAF files.")
    parser.add_argument("-s", "--specie", type=str, required=True, help="Species name used in the MAF file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Directory to save output files.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for parallel processing.")

    args = parser.parse_args()

    # Validate input paths
    if not os.path.isfile(args.bed12):
        logging.error(f"BED12 file not found: {args.bed12}")
        return
    if not os.path.isdir(args.dir):
        logging.error(f"Directory not found: {args.dir}")
        return
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        logging.info(f"Created output directory: {args.output}")

    try:
        bed12_entries = parse_bed12(args.bed12)

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            future_to_entry = {
                executor.submit(maf_parse, args.dir, entry, args.specie): entry for entry in bed12_entries
            }
            for future in as_completed(future_to_entry):
                entry = future_to_entry[future]
                try:
                    future.result()
                except Exception as e:
                    logging.error(f"Error processing entry {entry['name']}: {e}")

    except Exception as e:
        logging.error(f"Error during processing: {e}")

if __name__ == "__main__":
    main()
