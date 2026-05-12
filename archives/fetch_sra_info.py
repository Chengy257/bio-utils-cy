#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /home/chengyu/myscripts/fetch_sra_info.py
# Author: ChengYu
# Description: 
# Created Time: Sat 15 Nov 2025 02:49:55 PM CST
#########################################################################


import argparse
import csv
from Bio import Entrez
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import islice
import time

parser = argparse.ArgumentParser(description="Fetch SRA info and output SraRunTable style CSV")
parser.add_argument("-i", "--input", required=True, help="Input file with SRA IDs")
parser.add_argument("-o", "--output", default="SraRunTable.csv", help="Output CSV file")
parser.add_argument("-e", "--email", default="1067338303@qq.com", help="Email for NCBI Entrez API")
parser.add_argument("-b", "--batch", type=int, default=10, help="Batch size per request")
parser.add_argument("-t", "--threads", type=int, default=3, help="Number of parallel threads")
args = parser.parse_args()

Entrez.email = args.email

with open(args.input) as f:
    sra_ids = [line.strip() for line in f if line.strip()]

# ---------------------------
# Field names
# ---------------------------
fields = [
    "Run","BioProject","BioSample","Sample_Name","Organism","Tissue","Age","Cultivar",
    "Dev_stage","Ecotype","Consent","Geo_loc_country","Geo_loc_continent","Geo_loc",
    "Experiment","Library_Name","LibraryLayout","LibrarySource","LibrarySelection","Assay_Type",
    "Platform","Instrument","Center_Name","Spot_Count","Bases","AvgSpotLen",
    "DATASTORE_filetype","DATASTORE_provider","DATASTORE_region","ReleaseDate","CreateDate"
]

# ---------------------------
# Helper: chunk list
# ---------------------------
def chunks(lst, n):
    it = iter(lst)
    while True:
        c = list(islice(it, n))
        if not c:
            return
        yield c

# ---------------------------
# Parse one EXPERIMENT_PACKAGE
# ---------------------------
def parse_experiment_package(exp_pkg):
    rows = []
    # Study / BioProject
    study = exp_pkg.find("STUDY")
    bioproject = study.get("accession") if study is not None else ""

    # Sample
    sample = exp_pkg.find("SAMPLE")
    biosample = sample.findtext("IDENTIFIERS/BIO_SAMPLE","") if sample is not None else ""
    sample_name = sample.findtext("IDENTIFIERS/PRIMARY_ID","") if sample is not None else ""
    organism = sample.findtext("SAMPLE_NAME/SCIENTIFIC_NAME","") if sample is not None else ""

    # Sample attributes
    attr_map = {field: "" for field in ["Tissue","Age","Cultivar","Dev_stage","Ecotype","Consent",
                                        "Geo_loc_country","Geo_loc_continent","Geo_loc"]}
    if sample is not None:
        attrs = sample.find("SAMPLE_ATTRIBUTES")
        if attrs is not None:
            for a in attrs.findall("SAMPLE_ATTRIBUTE"):
                tag = a.findtext("TAG","").lower()
                val = a.findtext("VALUE","")
                if tag in ["tissue","age","cultivar","dev_stage","ecotype","consent",
                           "geo_loc_name_country","geo_loc_name_country_continent","geo_loc_name"]:
                    # 对应字段名
                    mapping = {
                        "tissue":"Tissue","age":"Age","cultivar":"Cultivar","dev_stage":"Dev_stage",
                        "ecotype":"Ecotype","consent":"Consent",
                        "geo_loc_name_country":"Geo_loc_country",
                        "geo_loc_name_country_continent":"Geo_loc_continent",
                        "geo_loc_name":"Geo_loc"
                    }
                    attr_map[mapping[tag]] = val

    # Experiment
    exp = exp_pkg.find("EXPERIMENT")
    experiment_acc = exp.get("accession","") if exp is not None else ""
    create_date = exp.get("submission_date","") if exp is not None else ""
    release_date = exp.get("published","") if exp is not None else ""

    # Library
    lib = exp.find("DESIGN/LIBRARY_DESCRIPTOR") if exp is not None else None
    lib_name = lib.findtext("LIBRARY_NAME","") if lib is not None else ""
    library_layout = "PAIRED" if lib.find("LIBRARY_LAYOUT/PAIRED") is not None else "SINGLE" if lib is not None else ""
    library_source = lib.findtext("LIBRARY_SOURCE","") if lib is not None else ""
    library_selection = lib.findtext("LIBRARY_SELECTION","") if lib is not None else ""
    assay_type = lib.findtext("LIBRARY_STRATEGY","") if lib is not None else ""

    # Platform / Instrument / Center
    plat = exp.find("PLATFORM/ILLUMINA") if exp is not None else None
    platform = plat.findtext("INSTRUMENT_MODEL","") if plat is not None else ""
    instrument = platform
    center_name = plat.findtext("CENTER_NAME","") if plat is not None else ""

    # Runs
    run_set = exp_pkg.findall(".//RUN")
    for run in run_set:
        row = {
            "Run": run.get("accession",""),
            "BioProject": bioproject,
            "BioSample": biosample,
            "Sample_Name": sample_name,
            "Organism": organism,
            **attr_map,
            "Experiment": experiment_acc,
            "Library_Name": lib_name,
            "LibraryLayout": library_layout,
            "LibrarySource": library_source,
            "LibrarySelection": library_selection,
            "Assay_Type": assay_type,
            "Platform": platform,
            "Instrument": instrument,
            "Center_Name": center_name,
            "Spot_Count": run.get("total_spots",""),
            "Bases": run.get("total_bases",""),
            "AvgSpotLen": run.get("avgSpotLen",""),
            "DATASTORE_filetype": run.findtext("DATASTORE/FILETYPE",""),
            "DATASTORE_provider": run.findtext("DATASTORE/PROVIDER",""),
            "DATASTORE_region": run.findtext("DATASTORE/REGION",""),
            "ReleaseDate": release_date,
            "CreateDate": create_date
        }
        rows.append(row)
    return rows

# ---------------------------
# Fetch batch
# ---------------------------
def fetch_batch(batch_ids):
    try:
        ids_str = ",".join(batch_ids)
        handle = Entrez.efetch(db="sra", id=ids_str, rettype="xml")
        xml_data = handle.read()
        handle.close()
        root = ET.fromstring(xml_data)
        all_rows = []
        for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            all_rows.extend(parse_experiment_package(pkg))
        time.sleep(0.34)  # NCBI rate limit
        return all_rows
    except Exception as e:
        print(f"Warning: batch {batch_ids} failed: {e}")
        return []

# ---------------------------
# Parallel fetch
# ---------------------------
all_rows = []
with ThreadPoolExecutor(max_workers=args.threads) as executor:
    futures = [executor.submit(fetch_batch, chunk) for chunk in chunks(sra_ids, args.batch)]
    for f in as_completed(futures):
        all_rows.extend(f.result())

# ---------------------------
# Write CSV
# ---------------------------
with open(args.output,"w",newline="",encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for row in all_rows:
        writer.writerow(row)

print(f"Finished. Output saved to {args.output}")
