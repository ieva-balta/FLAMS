#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Merge UniProt CSV PTMs with EBI CSV PTMs
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import argparse
import requests
import os

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

ORGANISM_CACHE_FILE = "organism_taxid_cache.csv"



def load_taxid_cache():
    if os.path.exists(ORGANISM_CACHE_FILE):
        return pd.read_csv(ORGANISM_CACHE_FILE).set_index("taxid")["name"].to_dict()
    return {}


def save_taxid_cache(cache_dict):
    df = pd.DataFrame([{"taxid": k, "name": v} for k, v in cache_dict.items()])
    df.to_csv(ORGANISM_CACHE_FILE, index=False)


taxid_cache = load_taxid_cache()


def taxid_to_name(taxid):
    """
    Convert NCBI TaxID to scientific name (cached)
    """
    taxid = str(taxid)

    if taxid in taxid_cache:
        return taxid_cache[taxid]

    url = f"https://rest.uniprot.org/taxonomy/{taxid}"
    try:
        r = requests.get(url, timeout=10)
        if not r.ok:
            taxid_cache[taxid] = taxid
            return taxid

        name = r.json().get("scientificName", taxid)
        taxid_cache[taxid] = name
        save_taxid_cache(taxid_cache)
        return name

    except Exception:
        taxid_cache[taxid] = taxid
        return taxid



def harmonize_uniprot(df):
    """
    Standardize UniProt CSV into unified schema
    """

    df = df.copy()

    df["Entry_type"] = "UniProt"
    df["Organism"] = df["Organism"].fillna("Unknown")

    # Rename UniProt fields into unified scheme
    df = df.rename(columns={
        "ECO_codes": "Evidence",
        "Sources": "Source_types",
        "Source_ids": "Source_ids",
        "Dataset": "Dataset"
    })

    df["Source_types"] = df["Source_types"].fillna("")
    df["Source_ids"] = df["Source_ids"].fillna("")
    df["Evidence"] = df["Evidence"].fillna("")
    df["Description"] = df["Description"].fillna("")

    standard_cols = [
    "Unique_ID", "Accession", "Position", "Length", "Protein", "Description",
    "Organism", "Evidence", "Source_types", "Source_ids", "Sequence",
    "Entry_type"
    ]


    return df[standard_cols]


def harmonize_ebi(df):
    """
    Standardize EBI CSV into unified schema
    """

    df = df.copy()

    # convert taxid to scientific name
    df["Organism"] = df["Organism"].astype(str).apply(taxid_to_name)

    df["Entry_type"] = "EBI"

    df = df.rename(columns={
        "ECO_codes": "Evidence",
        "Sources": "Source_types",
        "Source_ids": "Source_ids",
        "Description": "Description",
    })

    standard_cols = [
    "Unique_ID", "Accession", "Position", "Length", "Protein", "Description",
    "Organism", "Evidence", "Source_types", "Source_ids", "Sequence",
    "Entry_type"
    ]

    return df[standard_cols]


# def fix_source_order(source_types, source_ids):
#     """
#     Ensure Source_ids follows the SAME ORDER as Source_types,
#     and that PXD always precedes PubMed where appropriate.
#     """

#     # Split
#     dbs = [s for s in source_types.split(";") if s]
#     ids = [s for s in source_ids.split(";") if s]

#     # Classify
#     pxd_ids = [i for i in ids if i.startswith("PXD")]
#     pubmed_ids = [i for i in ids if i.isdigit()]

#     ordered_ids = []

#     for db in dbs:
#         if db == "PXD":
#             ordered_ids.extend(pxd_ids)
#         elif db == "PubMed":
#             ordered_ids.extend(pubmed_ids)

#     return ";".join(ordered_ids)


def merge_unique(column):
    """
    Deduplicate and merge unique items in a column
    """
    items = set()
    for v in column:
        if pd.isna(v):
            continue
        for x in str(v).split(";"):
            x = x.strip()
            if x:
                items.add(x)
    return ";".join(sorted(items))


# def deduplicate(df):

#     # Group rows by site (Accession + Position)
#     grouped = (
#         df.groupby(["Accession", "Position"], as_index=False)
#         .agg({
#             "Length": "first",        # assume same length
#             "Protein": "first",       # assume same protein
#             "Organism": "first",      # might differ but we keep first
#             "Description": merge_unique,     # merge all PTM types
#             "Evidence": merge_unique,        # merge evidence codes
#             "Source_types": merge_unique,    # merge: PXD;PubMed;...
#             "Source_ids": merge_unique,      # merge: PXDxxxxx;PubMed…
#             "Sequence": "first",             # identical for all rows
#         })
#     )

#     return grouped

    # Now enforce ordering of Source_ids to match Source_types
    # def reorder(row):
    #     row["Source_ids"] = fix_source_order(
    #         row["Source_types"], row["Source_ids"]
    #     )
    #     return row

    # grouped = grouped.apply(reorder, axis=1)

    # return grouped



def df_to_fasta(df, out_path):
    records = []

    for _, r in df.iterrows():
        seq = Seq(r["Sequence"])
        desc = (
            f"{r['Protein']}|{r['Description']}|{r['Organism']} "
            f"[{r['Source_types']}|{r['Source_ids']}|{r['Evidence']}]"
        )

        rec = SeqRecord(
            seq,
            id=f"{r['Accession']}|{r['Position']}|{r['Length']}|{r['Protein']}",
            description=desc
        )

        records.append(rec)

    SeqIO.write(records, out_path, "fasta")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniprot", required=True)
    parser.add_argument("--ebi", required=True)
    parser.add_argument("--out", default="merged_output")

    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    logging.info("Loading UniProt CSV…")
    uni_raw = pd.read_csv(args.uniprot)
    uni = harmonize_uniprot(uni_raw)

    logging.info("Loading EBI CSV…")
    ebi_raw = pd.read_csv(args.ebi)
    ebi = harmonize_ebi(ebi_raw)

    logging.info("Merging datasets")
    merged = pd.concat([uni, ebi], ignore_index=True)

    # logging.info("Concatenating datasets…")
    # combined = pd.concat([uni, ebi], ignore_index=True)

    # logging.info("Deduplicating merged dataset…")
    # merged = deduplicate(combined)


    out_csv = f"{args.out}/uniprot_ebi_merged.csv"
    out_fasta = f"{args.out}/uniprot_ebi_merged.fasta"

    merged.to_csv(out_csv, index=False)
    df_to_fasta(merged, out_fasta)

    logging.info(f"Saved merged CSV {out_csv}")
    logging.info(f"Saved merged FASTA {out_fasta}")


if __name__ == "__main__":
    main()
