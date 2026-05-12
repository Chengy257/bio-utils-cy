#!/usr/bin/env python3
#########################################################################
# File Name: batch_pride_fetch.py
# Author: ChengYu
# Description: 
# Created Time: Fri 21 Nov 2025 11:09:05 AM CST
#########################################################################


import requests
import sys

def fetch_pxd_info(pxd):
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{pxd}"
    try:
        r = requests.get(url, timeout=20)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        return pxd, "NA", f"ERROR: {e}"

    title = data.get("title", "NA")
    date = data.get("publicationDate") or data.get("submissionDate") or "NA"
    return pxd, date, title


def main():
    if len(sys.argv) != 3:
        print("Usage: python batch_pride_fetch.py pxd_list.txt output.tsv")
        sys.exit(1)

    pxd_file = sys.argv[1]
    out_file = sys.argv[2]

    # 读取 PXD 列表
    with open(pxd_file) as f:
        pxd_list = [line.strip() for line in f if line.strip()]

    with open(out_file, "w") as out:
        out.write("PXD_ID\tPublication_Date\tTitle\n")
        for pxd in pxd_list:
            pxd_id, date, title = fetch_pxd_info(pxd)
            out.write(f"{pxd_id}\t{date}\t{title}\n")
            print(f"Fetched: {pxd_id}")

    print(f"\nDone! Output saved to: {out_file}")


if __name__ == "__main__":
    main()
