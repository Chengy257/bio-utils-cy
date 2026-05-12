#!/usr/bin/env python3

import argparse
import requests
import pandas as pd
import time
import sys
import xml.etree.ElementTree as ET
from io import StringIO
import os


BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_EMAIL = os.getenv("NCBI_EMAIL", "your_email@lab.edu")
DEFAULT_API_KEY = os.getenv("NCBI_API_KEY")


def request_retry(url, params, retry=3):

    for i in range(retry):

        try:

            r = requests.get(url, params=params, timeout=60)

            r.raise_for_status()

            return r.text

        except Exception:

            if i < retry-1:

                time.sleep(3)

            else:

                raise


def fetch_runinfo(ids, email, api_key):

    params = {

        "db": "sra",
        "rettype": "runinfo",
        "retmode": "text",
        "id": ",".join(ids),
        "email": email

    }

    if api_key:
        params["api_key"] = api_key

    text = request_retry(
        f"{BASE_URL}/efetch.fcgi",
        params
    )

    return pd.read_csv(StringIO(text))


def fetch_esummary(ids, email, api_key):

    params = {

        "db": "sra",
        "id": ",".join(ids),
        "retmode": "xml",
        "email": email

    }

    if api_key:
        params["api_key"] = api_key

    xml = request_retry(
        f"{BASE_URL}/esummary.fcgi",
        params
    )

    return xml


def parse_xml_metadata(xml):

    root = ET.fromstring(xml)

    rows = []

    for doc in root.findall(".//DocSum"):

        data = {}

        for item in doc.findall("Item"):

            name = item.attrib.get("Name")

            if name == "ExpXml":

                try:

                    expxml = item.text

                    x = ET.fromstring(expxml)

                    # experiment info
                    exp = x.find(".//EXPERIMENT")

                    if exp is not None:

                        data["Experiment"] = exp.attrib.get("accession")

                        data["LibraryStrategy"] = exp.findtext(".//LIBRARY_STRATEGY")

                        data["LibrarySource"] = exp.findtext(".//LIBRARY_SOURCE")

                        data["LibrarySelection"] = exp.findtext(".//LIBRARY_SELECTION")

                        data["LibraryLayout"] = (
                            "PAIRED"
                            if exp.find(".//PAIRED") is not None
                            else "SINGLE"
                        )

                    # platform
                    platform = x.find(".//PLATFORM")

                    if platform is not None:

                        for child in platform:

                            data["Platform"] = child.tag

                            data["Model"] = child.text

                    # sample info
                    sample = x.find(".//SAMPLE")

                    if sample is not None:

                        data["SampleName"] = sample.attrib.get("alias")

                        data["SampleAccession"] = sample.attrib.get("accession")

                        organism = sample.find(".//SCIENTIFIC_NAME")

                        if organism is not None:

                            data["Organism"] = organism.text

                        # sample attributes
                        for attr in sample.findall(".//SAMPLE_ATTRIBUTE"):

                            tag = attr.findtext("TAG")

                            val = attr.findtext("VALUE")

                            if tag and val:

                                data[tag] = val

                    # study
                    study = x.find(".//STUDY")

                    if study is not None:

                        data["Study"] = study.attrib.get("accession")

                        title = study.findtext(".//STUDY_TITLE")

                        if title:
                            data["StudyTitle"] = title

                except:
                    pass

        rows.append(data)

    return pd.DataFrame(rows)


def read_ids(file):

    ids = []

    with open(file) as f:

        for line in f:

            line=line.strip()

            if line:
                ids.append(line)

    unique = list(dict.fromkeys(ids))

    if len(ids) != len(unique):

        print(f"Removed {len(ids)-len(unique)} duplicate IDs")

    return unique


def chunk(lst, n):

    for i in range(0, len(lst), n):

        yield lst[i:i+n]


def main():

    parser = argparse.ArgumentParser(
        description="Fetch comprehensive SRA metadata"
    )

    parser.add_argument("-i","--input",required=True)

    parser.add_argument("-o","--output",default="sra_metadata.tsv")

    parser.add_argument("--email",default=DEFAULT_EMAIL)

    parser.add_argument("--api-key",default=DEFAULT_API_KEY)

    parser.add_argument("--batch",type=int,default=200)

    args = parser.parse_args()

    print("NCBI configuration:")

    print("Email:",args.email)

    print("API key:", "provided" if args.api_key else "None")

    ids = read_ids(args.input)

    print("Total IDs:",len(ids))

    run_tables=[]
    xml_tables=[]

    for group in chunk(ids,args.batch):

        print("Querying",len(group),"IDs...")

        run_df = fetch_runinfo(
            group,
            args.email,
            args.api_key
        )

        xml = fetch_esummary(
            group,
            args.email,
            args.api_key
        )

        xml_df = parse_xml_metadata(xml)

        run_tables.append(run_df)

        xml_tables.append(xml_df)

        time.sleep(0.34)

    runinfo = pd.concat(run_tables)

    xmldata = pd.concat(xml_tables)

    meta = pd.concat(
        [runinfo.reset_index(drop=True),
         xmldata.reset_index(drop=True)],
        axis=1
    )

    meta.to_csv(
        args.output,
        sep="\t",
        index=False
    )

    print("\nSaved to:",args.output)

    print("Total runs:",len(meta))


if __name__=="__main__":
    main()