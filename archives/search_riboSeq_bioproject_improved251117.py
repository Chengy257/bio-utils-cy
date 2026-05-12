#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
search_riboseq_bioproject_recall_v2.py
--------------------------------------
Ribo-seq & Translatome & Polysome 项目的最大召回 BioProject 检索脚本

特性：
  ✓ 最强关键词词典（strict/medium/broad/fallback）
  ✓ 多层检索 + 合并去重
  ✓ 物种名 → taxid 自动识别
  ✓ BioProject → SRA 自动关联
  ✓ CSV 输出

依赖：
  pip install biopython
"""

import argparse
import csv
import time
import sys
from typing import List, Dict, Optional
from urllib.error import HTTPError, URLError
from Bio import Entrez
import xml.etree.ElementTree as ET

Entrez.email = "1067338303@qq.com"   # ⚠️ 修改为你的邮箱
SLEEP = 0.34


# ============================================================
# 关键词词典（最强版）
# ============================================================

STRICT_KEYWORDS = [
    "ribosome profiling",
    "ribo-seq", "riboseq", "ribo seq",
    "ribosome footprinting", "ribosome footprints",
    "ribosome protected fragments", "ribosome-protected fragments",
    "RPF", "RPF-seq"
]

MEDIUM_KEYWORDS = [
    # translatome
    "translatome", "translatomics",
    "translatomic profiling", "translatome sequencing",
    "translation-level profiling",
    "actively translated mRNA",
    "actively translating",

    # polysome profiling
    "polysome profiling", "polysome-seq", "polysome sequencing",
    "polyribosome profiling",
    "polysome fractionation",
    "ribosome-associated mRNA",
    "ribosome bound mRNA",
    "ribosome-associated transcriptome",

    # translation efficiency / TE
    "TE profiling",
    "translation-level analysis",

    # translation dynamics
    "ribosome occupancy",
    "ribosome release",
    "ribosome pausing",
    "ribosome stalling",
    "translation dynamics"
]

# 宽泛匹配（补充召回）
BROAD_PATTERNS = [
    '(ribosome[All Fields] AND (profiling[All Fields] OR footprinting[All Fields]))',
    '(translation[All Fields] AND (profiling[All Fields] OR footprint[All Fields]))'
]

# fallback（最宽泛）
# FALLBACK_PATTERN = '(ribosome[All Fields] AND sequencing[All Fields])'


# ============================================================
# taxonomy 查询
# ============================================================
def species_to_taxid(species: str) -> Optional[int]:
    try:
        h = Entrez.esearch(db="taxonomy", term=f'"{species}"[Scientific Name]')
        data = Entrez.read(h)
        h.close()
        time.sleep(SLEEP)
        ids = data.get("IdList", [])
        if ids:
            return int(ids[0])
    except Exception as e:
        print(f"[WARN] taxonomy lookup failed for {species}: {e}", file=sys.stderr)
    return None


# ============================================================
# 多层查询构建
# ============================================================
def build_queries(taxid: Optional[int]) -> List[str]:
    suffix = f" AND txid{taxid}[Organism]" if taxid else ""

    # Layer 1: STRICT
    strict_query = " OR ".join([f'"{k}"[All Fields]' for k in STRICT_KEYWORDS])

    # Layer 2: MEDIUM
    medium_query = " OR ".join([f'"{k}"[All Fields]' for k in MEDIUM_KEYWORDS])

    # Layer 3: BROAD
    broad_query = " OR ".join(BROAD_PATTERNS)

    # Layer 4: FALLBACK
    # fallback_query = FALLBACK_PATTERN

    return [
        f"({strict_query}){suffix}",
        f"({medium_query}){suffix}",
        f"({broad_query}){suffix}",
        # f"({fallback_query}){suffix}"
    ]


# ============================================================
# esearch（全量检索）
# ============================================================
def esearch_all(query: str, retmax=500) -> List[str]:
    ids = []
    retstart = 0
    print(f"[QUERY] {query}", file=sys.stderr)

    while True:
        try:
            h = Entrez.esearch(
                db="bioproject",
                term=query,
                retmax=retmax,
                retstart=retstart,
                usehistory="n"
            )
            data = Entrez.read(h)
            h.close()
            time.sleep(SLEEP)
        except Exception as e:
            print(f"[ERROR] esearch: {e}", file=sys.stderr)
            break

        batch = data.get("IdList", [])
        if not batch:
            break

        ids.extend(batch)
        if len(batch) < retmax:
            break

        retstart += retmax

    return ids


# ============================================================
# efetch BioProject details
# ============================================================

# def fetch_bioproject(ids: List[str]) -> List[Dict]:
#     if not ids:
#         return []
#     out = []
#     chunk = 50

#     for i in range(0, len(ids), chunk):
#         sub = ids[i:i + chunk]
#         try:
#             h = Entrez.efetch(db="bioproject", id=",".join(sub), rettype="xml")
#             xml = h.read()
#             h.close()
#             time.sleep(SLEEP)
#         except Exception as e:
#             print(f"[WARN] efetch error: {e}", file=sys.stderr)
#             continue

#         try:
#             root = ET.fromstring(xml)
#         except Exception as e:
#             print(f"[WARN] XML parse error: {e}", file=sys.stderr)
#             continue

#         for p in root.findall(".//Project"):
#             d = {}
#             arch = p.find("ProjectID/ArchiveID")
#             d["Accession"] = arch.attrib.get("accession", "") if arch is not None else ""

#             pid = p.find("ProjectID/ProjectID")
#             d["BioProjectID"] = pid.text if pid is not None else ""

#             title = p.find("Description/Title")
#             d["Title"] = title.text if title is not None else ""

#             org = p.find("Organism/OrganismName")
#             d["Organism"] = org.text if org is not None else ""

#             desc = p.find("Description/ProjectDescription")
#             d["Description"] = desc.text if desc is not None else ""

#             rel = p.find("ProjectStatus/ReleaseDate")
#             d["ReleaseDate"] = rel.text if rel is not None else ""

#             out.append(d)

#     return out

def fetch_bioproject(ids: List[str]) -> List[Dict]:
    """
    强健版解析：自动适配 NCBI / ENA / DDBJ 三种 BioProject XML 结构
    """
    if not ids:
        return []
    out = []
    chunk = 50

    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        try:
            h = Entrez.efetch(db="bioproject", id=",".join(sub), rettype="xml")
            xml = h.read()
            h.close()
            time.sleep(SLEEP)
        except Exception as e:
            print(f"[WARN] efetch error: {e}", file=sys.stderr)
            continue

        try:
            root = ET.fromstring(xml)
        except Exception as e:
            print(f"[WARN] XML parse error: {e}", file=sys.stderr)
            continue

        # 遍历所有 Project
        for p in root.findall(".//Project"):
            d = {}

            # === Accession（唯一绝对存在）===
            arch = p.find(".//ArchiveID")
            d["Accession"] = arch.attrib.get("accession", "") if arch is not None else ""

            # === BioProject numeric ID（一部分项目无该字段）===
            pid = p.find(".//ProjectID")
            if pid is not None and pid.text:
                d["BioProjectID"] = pid.text.strip()
            else:
                d["BioProjectID"] = d["Accession"]  # fallback

            # === Title（多种位置）===
            title_paths = [
                ".//Title",
                ".//ProjectDescr/Title",
                ".//Description/Title",
                ".//Name",  # 某些 ENA/DDBJ 用这个
            ]
            title = None
            for path in title_paths:
                node = p.find(path)
                if node is not None and node.text:
                    title = node.text.strip()
                    break
            d["Title"] = title if title else ""

            # === Organism ===
            organism_paths = [
                ".//Organism/OrganismName",
                ".//Organism/Name",
                ".//Organism"
            ]
            organism = None
            for path in organism_paths:
                node = p.find(path)
                if node is not None and node.text:
                    organism = node.text.strip()
                    break
            d["Organism"] = organism if organism else ""

            # === Description ===
            desc_paths = [
                ".//ProjectDescr/Description",
                ".//ProjectDescription",
                ".//Description/ProjectDescription",
                ".//Description",
            ]
            desc = None
            for path in desc_paths:
                node = p.find(path)
                if node is not None and node.text:
                    desc = node.text.strip()
                    break
            d["Description"] = desc if desc else ""

            # === Release Date ===
            date_paths = [
                ".//ReleaseDate",
                ".//Submission/ReleaseDate",
                ".//ProjectReleaseDate"
            ]
            rel = None
            for path in date_paths:
                node = p.find(path)
                if node is not None and node.text:
                    rel = node.text.strip()
                    break
            d["ReleaseDate"] = rel if rel else ""

            out.append(d)

    return out




# ============================================================
# BioProject → SRA
# ============================================================
def link_sra(bid: str) -> List[str]:
    out = []
    try:
        h = Entrez.elink(dbfrom="bioproject", db="sra", id=bid)
        data = Entrez.read(h)
        h.close()
        time.sleep(SLEEP)
        for block in data:
            for db in block.get("LinkSetDb", []):
                if db.get("DbTo") == "sra":
                    for link in db.get("Link", []):
                        out.append(link.get("Id"))
    except Exception:
        pass
    return out


# ============================================================
# CSV 保存
# ============================================================
def save_csv(records: List[Dict], outfile: str):
    fields = [
        "QuerySpecies", "BioProjectID", "Accession",
        "Title", "Organism", "ReleaseDate", "Description", "LinkedSRA"
    ]
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        for r in records:
            r["LinkedSRA"] = ";".join(r.get("LinkedSRA", []))
            w.writerow(r)


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Ribo-seq 强召回 BioProject 检索")
    parser.add_argument("-s", "--species", nargs="+", help="species name or taxid", default=[])
    parser.add_argument("-o", "--output", default="ribo_bioproject_recall_v2.csv")
    args = parser.parse_args()

    species_list = args.species if args.species else [None]
    final_records = []

    for sp in species_list:
        # 获取 taxid
        taxid = None
        label = ""

        if sp:
            try:
                taxid = int(sp)
                label = f"taxid:{taxid}"
            except ValueError:
                tx = species_to_taxid(sp)
                taxid = tx
                label = f"{sp}(txid={taxid})"
        else:
            label = "ALL_SPECIES"

        print(f"\n==== 物种：{label} ====", file=sys.stderr)

        # 多层查询
        queries = build_queries(taxid)

        id_set = set()
        for q in queries:
            ids = esearch_all(q)
            id_set.update(ids)

        print(f"[INFO] Unique BioProject IDs: {len(id_set)}", file=sys.stderr)

        # fetch
        details = fetch_bioproject(list(id_set))

        # link SRA
        for d in details:
            bid = d.get("BioProjectID", "")
            d["LinkedSRA"] = link_sra(bid)
            d["QuerySpecies"] = label
            final_records.append(d)

    # 保存
    save_csv(final_records, args.output)
    print(f"[DONE] Saved {len(final_records)} records → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
