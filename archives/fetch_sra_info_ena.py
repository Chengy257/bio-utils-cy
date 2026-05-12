#########################################################################
# File Name: /home/chengyu/myscripts/fetch_sra_info_ena.py
# Author: ChengYu
# Description: 
# Created Time: Sat 15 Nov 2025 03:20:12 PM CST
#########################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import requests
import pandas as pd
import io

parser = argparse.ArgumentParser(description="Fetch SRA info from ENA filereport API")
parser.add_argument("-i", "--input", required=True, help="Input file with SRA IDs")
parser.add_argument("-o", "--output", default="SraRunInfo.tsv", help="Output TSV file")
parser.add_argument("-f", "--fields", default="study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_name,nominal_length,library_selection,first_public,last_updated,experiment_title,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_ftp,sra_ftp,sample_alias",
                    help="Comma-separated list of fields to fetch")
args = parser.parse_args()

# Read SRA IDs
with open(args.input) as f:
    sra_ids = [line.strip() for line in f if line.strip()]

# Prepare list to collect dataframes
dfs = []

for sra_id in sra_ids:
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={sra_id}&result=read_run&fields={args.fields}&format=tsv&download=true&limit=0"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        df = pd.read_csv(io.StringIO(r.text), sep="\t")  # <--- 修改这里
        dfs.append(df)
        print(f"Fetched {sra_id}, {len(df)} rows")
    except Exception as e:
        print(f"Warning: failed to fetch {sra_id}: {e}")

# Merge all results
if dfs:
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(args.output, sep="\t", index=False)
    print(f"Finished. Output saved to {args.output}")
else:
    print("No data fetched.")
