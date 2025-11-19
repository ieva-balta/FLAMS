#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Downloader of large-scale PTMs from EBI Proteomics API using batch retrieval.
Uses UniProt existence:1 accessions.

Key features:
 - Streams results directly into CSV/FASTA (no large in-memory data structures)
 - Merges duplicate PTM sites per accession
 - Extracts ECO codes + PubMed IDs + Dataset IDs + PTM source (PRIDE / PTMeXchange)
 - Automatically resumes using processed.txt
 - Checkpoints every 5000 proteins
 - **BATCHES** up to 100 accessions per API request
"""

import argparse
import logging
import os
import requests
import pandas as pd
import re
import random
import string
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
# Note: Changing endpoint to the Proteomics search endpoint for batch query
EBI_PTM_SEARCH_URL = "https://www.ebi.ac.uk/proteins/api/proteomics/ptm"
HEADERS = {"Accept": "application/json"}
UNIPROT_PAGE_SIZE = 500
BATCH_SIZE = 100


def init_logging(outdir):
    log_file = os.path.join(outdir, "ebi_stream.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging to {log_file}")


# def is_valid_accession(acc):
#     """
#     Ensures the string is a valid UniProt accession
#     """
#     if not isinstance(acc, str):
#         return False
#     acc = acc.strip()

#     # Swiss-Prot: 1 letter (O,P,Q) + digit + 3 alphanumerics + digit
#     if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$", acc):
#         return True

#     # General UniProt format: 6 alphanumerics
#     if re.match(r"^[A-NR-Z0-9]{6}$", acc):
#         return True

#     return False


def fetch_uniprot_accessions():
    # ... (No change to this function, it fetches all accessions)
    # This is the url
    #https://rest.uniprot.org/uniprotkb/search?query=existence:1&fields=accession&includeIsoform=true&format=json&size=500
    url = f"{UNIPROT_URL}?query=existence:1&fields=accession&includeIsoform=true&format=json&size={UNIPROT_PAGE_SIZE}"
    accessions = []

    logging.info("Fetching UniProt accessions (existence:1)")
    while url:
        r = requests.get(url)
        r.raise_for_status()
        data = r.json()

        for entry in data.get("results", []):
            acc = entry.get("primaryAccession")
            if acc:
                accessions.append(acc)
        
        # lines = r.text.splitlines()
        # if lines and lines[0].lower() == "accession":
        #     lines = lines[1:]

        # accessions.extend(lines)

        next_link = r.links.get("next")
        url = next_link["url"] if next_link else None

        logging.info(f"Fetched {len(accessions)} accessions so far...")

    logging.info(f"Total existence:1 accessions: {len(accessions)}")
    return accessions


def batch_accessions(accessions, size):
    """Yields batches of accessions from a list."""
    for i in range(0, len(accessions), size):
        yield accessions[i:i + size]


# get the scientific name from taxid
TAX_CACHE = {}

def taxid_to_name(taxid):
    """
    Convert NCBI TaxID to scientific name using UniProt REST API
    """
    # ... (No change to this function)
    taxid = str(taxid)

    if taxid in TAX_CACHE:
        return TAX_CACHE[taxid]

    url = f"https://rest.uniprot.org/taxonomy/{taxid}"
    try:
        r = requests.get(url, timeout=10)
        if not r.ok:
            TAX_CACHE[taxid] = taxid
            return taxid
        name = r.json().get("scientificName", taxid)
        TAX_CACHE[taxid] = name.replace(" ", "__")
        return TAX_CACHE[taxid]
    except Exception:
        TAX_CACHE[taxid] = taxid
        return taxid



def random_tag(n=4):
    """
    Return a random uppercase letter combination as a tag
    """
    return ''.join(random.choices(string.ascii_uppercase, k=n))


# parse ptm data for a single accession (now handles a full batch response)
def process_ptm_data(data):
    """
    Processes the JSON response for a single protein entry from the EBI PTM API.
    Returns a list of rows for CSV/FASTA output.
    """
    if "accession" not in data:
        return []

    acc = data.get("accession")
    seq = data.get("sequence", "")
    protein = data.get("entryName", "")
    length = len(seq)
    organism = data.get("taxid", "")

    merged = {}  # key: (pos, desc)

    for feature in data.get("features", []):
        if feature.get("type") != "PROTEOMICS_PTM":
            continue

        begin = feature.get("begin")
        if begin is None:
            continue

        begin = int(begin)

        # ECO codes
        eco_codes = {
            ev.get("code")
            for ev in feature.get("evidences", [])
            if ev.get("code")
        }

        for ptm in feature.get("ptms", []):
            desc = ptm.get("name", "")
            rel_pos = ptm.get("position")
            if rel_pos is None:
                continue

            pos = begin + rel_pos - 1

            # source (PRIDE/PTMeXchange)
            sources = set(ptm.get("sources", []))

            # dataset id (PXD)
            dataset_ids = {
                d.get("id", "")
                for d in ptm.get("dbReferences", [])
                if d.get("id")
            }

            # pubmed id
            pubmed_ids = {
                d.get("properties", {}).get("Pubmed ID", "")
                for d in ptm.get("dbReferences", [])
                if d.get("properties", {}).get("Pubmed ID")
            }

            key = (pos, desc)
            if key not in merged:
                merged[key] = {
                    "Accession": acc,
                    "Position": pos,
                    "Length": length,
                    "Protein": protein,
                    "Description": desc,
                    "Organism": organism,
                    "Sources": sources,
                    "Dataset": dataset_ids,
                    "ECO": eco_codes,
                    "PubMed": pubmed_ids,
                    "Sequence": seq,
                }
            else:
                e = merged[key]
                e["Sources"].update(sources)
                e["Dataset"].update(dataset_ids)
                e["ECO"].update(eco_codes)
                e["PubMed"].update(pubmed_ids)

    # convert merged dictionary to list of rows
    rows = []
    for v in merged.values():
        rows.append({
            # "Unique_ID": f"{v['Accession']}_{v['Position']}",
            "Unique_ID": f"{v['Accession']}_{random_tag()}",
            "Accession": v["Accession"],
            "Position": v["Position"],
            "Length": v["Length"],

            "Protein": v["Protein"].replace(" ", "__"),
            "Feature_type": "Large_scale_PTM",
            "Description": v["Description"].replace(" ", "__"),
            "Organism": taxid_to_name(v["Organism"]),

            # database (PRIDE/PTMeXchange)
            "Database": ";".join(sorted(v["Sources"])),

            # sources (pdb, pubmed)
            "Sources": ";".join(sorted(
                ["PubMed" for _ in v["PubMed"]] +
                ["PXD" for _ in v["Dataset"]]
            )),

            # IDs of the above (PubMed IDs + PXDs)
            "Source_ids": ";".join(sorted(
                list(v["Dataset"]) + list(v["PubMed"])
            )),

            # ECO codes
            "ECO_codes": ";".join(sorted(v["ECO"])),

            "Sequence": v["Sequence"],
        })

    return rows


def fetch_ptms_for_batch(batch_accs):
    """
    Fetches PTMs for a batch of accessions from EBI Proteomics API.
    Handles pagination via 'next' links in response headers
    Returns a tuple of (list_of_PTM_rows, list_of_processed_accessions).
    """
    acc_list = ",".join(batch_accs)
    
    url = f"{EBI_PTM_SEARCH_URL}?accession={acc_list}&size=100" 
    
    # all_rows = []
    # processed_accs_in_response = set()
    
    # while url:
    #     r = requests.get(url, headers=HEADERS)
    #     if not r.ok:
    #         logging.error(f"API error: {r.status_code} - {r.text}")
    #         raise ValueError(f"API request failed for batch starting with {batch_accs[0]}")
    #     try:
    #         data = r.json()
    #     except ValueError:
    #         logging.error(f"Invalid JSON response for {url}")
    #         data = []
        
    #     for entry in data:
    #         rows = process_ptm_data(entry)
    #         all_rows.extend(rows)

    #         if "accession" in entry:
    #             processed_accs_in_response.add(entry["accession"])
        
    #     next_link = r.links.get("next")
    #     url = next_link["url"] if next_link else None
    #     logging.info(f"Processed page for batch; next URL: {url or 'None'}")

    r = requests.get(url, headers = HEADERS)
    if not r.ok:
        r.raise_for_status() # Raise error for a bad request
        return [], []

    # batch endpoint returns a list of protein entries
    data = r.json()
    all_rows = []
    processed_accs_in_response = set()
    
    # go through each entry in the batch
    for entry in data:
        rows = process_ptm_data(entry)
        all_rows.extend(rows)
        # mark accession as processed whether PTMs were found or not
        if "accession" in entry:
            processed_accs_in_response.add(entry["accession"])

    return all_rows, list(processed_accs_in_response)


# write to csv
def append_rows_to_csv(rows, csv_path):
    header_exists = os.path.isfile(csv_path)

    with open(csv_path, "a") as f:
        if not header_exists:
            f.write(
                "Unique_ID,Accession,Position,Description,Length,Database,Protein,"
                "Feature_type,Organism,ECO_codes,Sources,Source_ids,Sequence\n"
            )


        for r in rows:
            f.write(
                f"{r['Unique_ID']},{r['Accession']},{r['Position']},"
                f"{r['Description']},{r['Length']},{r['Database']},{r['Protein']},"
                f"{r['Feature_type']},{r['Organism']},{r['ECO_codes']},"
                f"{r['Sources']},{r['Source_ids']},{r['Sequence']}\n"
            )


def deduplicate(csv_path):
    if not os.path.isfile(csv_path):
        logging.warning("No CSV found for deduplication")
        return

    df = pd.read_csv(csv_path)

    grouped = (
        df.groupby(["Accession", "Position", "Description"], as_index=False)
        .agg({
            "Unique_ID": "first",
            "Length": "first",
            "Database": "first",
            "Protein": "first",
            "Feature_type": "first",
            "Organism": "first",
            "ECO_codes": lambda x: ";".join(sorted({i for s in x for i in str(s).split(";") if i})),
            "Sources":    lambda x: ";".join(sorted({i for s in x for i in str(s).split(";") if i})),
            "Source_ids": lambda x: ";".join(sorted({i for s in x for i in str(s).split(";") if i})),
            "Sequence": "first",
        })
    )

    correct_order = [
        "Unique_ID","Accession","Position","Description","Length","Database",
        "Protein","Feature_type","Organism","ECO_codes","Sources","Source_ids","Sequence"
    ]
    grouped = grouped[correct_order]

    grouped.to_csv(csv_path, index=False)
    logging.info(f"CSV overwritten with deduplicated results: {csv_path}")
    logging.info(f"After deduplication: {len(grouped)} rows")

    # write back to CSV
    # grouped.to_csv(csv_path, index=False)
    # logging.info(f"CSV overwritten with deduplicated results: {csv_path}")


def regenerate_fasta_from_csv(csv_path, fasta_path):
    df = pd.read_csv(csv_path)
    records = []

    for _, r in df.iterrows():
        seq = Seq(r["Sequence"])

        desc = (
            f"{r['Protein']}|{r['Description']}|{r['Organism']} "
            f"[{r['ECO_codes']}|{r['Sources']}|{r['Source_ids']}]"
        )

        rec = SeqRecord(
            seq,
            id=f"{r['Accession']}|{r['Position']}|{r['Length']}|{r['Database']}",
            description=desc
        )
        records.append(rec)

    SeqIO.write(records, fasta_path, "fasta")

    logging.info(f"Deduplicated FASTA regenerated: {fasta_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id", help="Run for single accession only (will be run as a batch of 1)")
    parser.add_argument("--out", default="ebi_output")
    parser.add_argument("--threads", type=int, default=40, help="Number of concurrent batches (max 100 accessions each)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    init_logging(args.out)

    csv_path = f"{args.out}/EBI_large_scale_PTMs.csv"
    fasta_path = f"{args.out}/EBI_large_scale_PTMs.fasta"
    processed_path = f"{args.out}/processed.txt"
    # no_ptm_path = f"{args.out}/no_ptm.accessions.txt"

    processed = set()
    if os.path.isfile(processed_path):
        with open(processed_path) as f:
            processed = {line.strip() for line in f}

    if args.id:
        accessions = [args.id]
        logging.info(f"SINGLE MODE: {args.id}")
    else:
        # accessions = fetch_uniprot_accessions()
        # raw_accs = fetch_uniprot_accessions()

        # Filter out invalid strings like "Entry", "Accession", empty, etc.
        accessions = fetch_uniprot_accessions()
        logging.info(f"Filtered to {len(accessions)} valid UniProt accessions.")


    to_process = [acc for acc in accessions if acc not in processed]
    logging.info(f"Total remaining accessions: {len(to_process)}")

    # Create batches of accessions
    batches = list(batch_accessions(to_process, BATCH_SIZE))
    logging.info(f"Total batches to process: {len(batches)}")


    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        # Submit each batch to the thread pool
        futures = {ex.submit(fetch_ptms_for_batch, batch): batch for batch in batches}

        new_batch_processed = 0


        for f in as_completed(futures):
            current_batch = futures[f]
            new_batch_processed += len(current_batch)
            try:
                rows, processed_accs = f.result()
                
                if rows:
                    append_rows_to_csv(rows, csv_path)
                    
                # Log success for the batch
                logging.info(
                    f"Batch processed ({new_batch_processed} accs): "
                    f"{len(rows)} PTM sites found."
                )

                # Update processed.txt with all accessions in the batch
                with open(processed_path, "a") as p:
                    for acc in current_batch:
                        p.write(acc + "\n")
                
                # Accessions without PTMs (not in API response)
                # no_ptm = set(current_batch) - set(processed_accs)
                # if no_ptm:
                #     with open(no_ptm_path, "a") as np:
                #         for acc in sorted(no_ptm):
                #             np.write(acc + "\n")
                #     logging.info(f"Added {len(no_ptm)} accessions without PTMs")

            except Exception as e:
                # Log an error for the batch but do NOT mark as processed
                # so it can be retried on the next run
                logging.error(
                    f"ERROR processing batch of {len(current_batch)} accessions "
                    f"starting with {current_batch[0]}: {e}"
                )


    logging.info("All PTMs fetched. Starting global deduplication...")
    deduplicate(csv_path)
    logging.info("CSV deduplicated. Generating FASTA")
    regenerate_fasta_from_csv(csv_path, fasta_path)
    logging.info("DONE.")


if __name__ == "__main__":
    main()