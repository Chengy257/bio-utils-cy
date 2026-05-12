#!/usr/bin/python
#########################################################################
# File Name: /home/chengyu/myscripts/calculate_average_pi.py
# Author: ChengYu
# Description: 
# Created Time: Sat 11 Jan 2025 05:27:16 PM CST
#########################################################################
#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor
import logging
import re

def setup_logger():
    """Set up logger for logging messages."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler("process.log", mode="w"),
            logging.StreamHandler()
        ]
    )

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Calculate SNP count and nucleotide diversity (Π) for genomic intervals.")
    parser.add_argument("-b", "--bed", required=True, help="Input BED6 file with genomic intervals.")
    parser.add_argument("-p", "--plink", required=True, help="PLINK prefix for SNP data (.bed, .bim, .fam).")
    parser.add_argument("-o", "--output", required=True, help="Output file to save results.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for parallel processing.")
    return parser.parse_args()

def sanitize_filename(name):
    """
    Sanitize the name to create a safe filename.
    Replace invalid characters with '_'.
    """
    return re.sub(r'[^\w.-]', '_', name)

def extract_snps_and_calculate_pi(plink_prefix, bed_row, output_dir):
    """Extract SNPs for a single BED interval and calculate site-level Pi."""
    chrom, start, end, name = bed_row

    # Clean name for file safety
    sanitized_name = sanitize_filename(name)
    logging.info(f"Processing interval: {name} ({chrom}:{start}-{end}), sanitized as: {sanitized_name}")

    # Define file paths
    range_file = os.path.join(output_dir, f"{sanitized_name}_range.txt")
    vcf_file = os.path.join(output_dir, f"{sanitized_name}.vcf")
    pi_file = os.path.join(output_dir, f"{sanitized_name}.site.pi")

    # Write range to file
    with open(range_file, "w") as f:
        f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    
    # Run PLINK to extract SNPs
    # plink_cmd = [
    #     "/home/chengyu/soft/miniconda3/envs/snp/bin/plink",
    #     "--bfile", plink_prefix,
    #     "--threads", "40",
    #     "--memory", "400000",
    #     "--extract 'range'", range_file,
    #     "--recode", "vcf-iid",
    #     "--out", vcf_file.replace(".vcf", "")
    # ]
    # print(" ".join(plink_cmd))
    # subprocess.run(plink_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Run vcftools to calculate site-level Pi
    vcftools_cmd = [
        "/home/chengyu/soft/miniconda3/envs/snp/bin/vcftools",
        "--vcf", vcf_file,
        "--site-pi",
        "--out", pi_file.replace(".site.pi", "")
    ]
    subprocess.run(vcftools_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Parse site-level Pi and calculate average
    pi_value = 0.0
    snp_count = 0
    try:
        with open(f"{pi_file}.site.pi", "r") as f:
            lines = f.readlines()[1:]  # Skip header
            if lines:
                values = [float(line.strip().split()[-1]) for line in lines]
                pi_value = sum(values) / len(values)  # Average Pi
                snp_count = len(values)
                logging.info(f"Interval {name}: SNP Count = {snp_count}, Average Pi = {pi_value:.6f}")
            else:
                logging.warning(f"Interval {name}: No SNPs found.")
    except FileNotFoundError:
        logging.warning(f"Interval {name}: No SNPs found (site-pi file missing).")
    
    return name, snp_count, pi_value

def process_bed_file(bed_file, plink_prefix, output_file, threads):
    """Process BED file and compute SNP count and average Pi values."""
    # Load BED file
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    results = []
    output_dir = "temp_output"
    os.makedirs(output_dir, exist_ok=True)

    # Run multithreading
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(extract_snps_and_calculate_pi, plink_prefix, row, output_dir)
            for row in bed_df.itertuples(index=False, name=None)
        ]
        for future in futures:
            results.append(future.result())

    # Save results to file
    results_df = pd.DataFrame(results, columns=["name", "snp_count", "pi_value"])
    results_df.to_csv(output_file, sep="\t", index=False)
    logging.info(f"Results saved to {output_file}")

if __name__ == "__main__":
    setup_logger()
    args = parse_arguments()
    process_bed_file(args.bed, args.plink, args.output, args.threads)
