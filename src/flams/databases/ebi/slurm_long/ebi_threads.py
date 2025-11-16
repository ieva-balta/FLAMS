#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import json
import pandas as pd
import csv
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import requests
import logging

import time
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


#valid eco codes
valid_ECO_codes = [
    "ECO:0000269", #(experimental evidence used in manual assertion)​
    "ECO:0000314", #(direct assay evidence used in manual assertion)​
    "ECO:0007744", #(combinatorial computational and experimental evidence used in manual assertion)​
    "ECO:0007829", #(combinatorial computational and experimental evidence used in automatic assertion)​
]


# columns for the csv of info 
columns = ["Accession", "Position", "Length", "Entry_type",
                            "Protein", "Feature_type", "Description", "Organism",
                            "ECO_codes", "Sources", # for all
                            "LSS_Database", "LSS_IDs", "LSS_Confidence_scores"# Only for large scale study data
                            ] 

# to store sequence separately
seq_columns = ["Accession", "Sequence"]




def process_entry(entry, headers):
    """Process a single entry from the EBI database and extract relevant information."""
    start_time = time.time()

    df_rows = []
    seq_rows = []

    accession = entry.get("primaryAccession")
    sequence = entry.get("sequence", {}).get("value", "")

    # stores sequence separately
    seq_rows.append([accession, sequence])

    name = entry.get("genes", [{}])[0].get("geneName", {}).get("value", "")
    organism = entry.get("organism", {}).get("scientificName", "")

    # storing database
    entry_type = entry.get("entryType", "")
    if "Swiss-Prot" in entry_type:
        entry_type = "Swiss-Prot"
    if "TrEMBL" in entry_type:
        entry_type = "TrEMBL"

    # removes white space
    entry_type = entry_type.replace(" ", "__")
    protein_name = f"{name}".replace(" ","__")
    organism_name = f"{organism}".replace(" ","__")
    length = len(sequence)
                
    # finds PTM sites from REST

    # to deduplicate the PTMs from REST
    ptms_rest = {}

    # finds relevant feature types
    for feature in entry.get("features", []):
        feature_type = feature.get("type")
        if feature.get("type") not in ["Modified residue",
                                # "Chain",
                                "Modified residue (large scale data)", 
                                "Lipidation", 
                                "Glycosylation",
                                "Disulfide bond",
                                "Cross-link"]:
            continue

        # PTM description
        desc = feature.get("description", "")

        #special case  - disulfide bonds
        if feature_type == "Disulfide bond":
            desc = "Disulfide bond"

        # remove the notes from the description
        if ";" in desc:
            desc = desc.split(";")[0].strip()

        #special case - "removed"
        if desc.lower() == "removed":
            continue

        #special case - microbial infection
        if "microbial infection" in desc.lower():
            desc = desc.replace("(Microbial infection) ", "")

        # gets positions for start and end (matters for crosslinks, bond formation, to record both amino acids)
        pos_start = feature.get("location", {}).get("start", {}).get("value", "N/A")
        pos_end = feature.get("location", {}).get("end", {}).get("value", "N/A")

        # gets evidence ECO|sources
        ECOs = set()
        sources = set()
        for ev in feature.get("evidences", []):
            eco = ev.get("evidenceCode", "")
            # skips if the evidnece is from the valid eco codes
            if eco not in valid_ECO_codes:
                continue
            source_site = ev.get("source", "")
            source_id = str(ev.get("id", ""))
            source = f"{source_site}:{source_id}" 
            ECOs.add(eco)
            sources.add(source)
                # ECO_str = ";".join(ECOs) 
                # # skips if there are no ECOs
                # if not ECOs:
                #     continue 
                # source_str = ";".join(sources) if sources else ""

        # removes white space
        feature_type_name = f"{feature_type}".replace(" ", "__")
        desc = f"{desc}".replace(" ","__")
                    
        # deduplicates based on position + description
        ptm_id = f"{pos_start}|{desc}"
        ptm_desc = [ECOs, sources, feature_type_name]

        if ptm_id not in ptms_rest:
            ptms_rest[ptm_id] = ptm_desc
        else:
            ptms_rest[ptm_id][0].update(ECOs)
            ptms_rest[ptm_id][1].update(sources)

        if pos_start != pos_end:
            ptm_id = f"{pos_end}|{desc}"
            if ptm_id not in ptms_rest:
                ptms_rest[ptm_id] = ptm_desc
            else:
                ptms_rest[ptm_id][0].update(ECOs)
                ptms_rest[ptm_id][1].update(sources)

    # appends to dataframe
    for key, item in ptms_rest.items():
        ecos_set = item[0]
        if not ecos_set:
            continue
        df_rows.append([accession, key.split("|")[0], length, entry_type,
                            protein_name, item[2], key.split("|")[1], organism_name,
                            ";".join(ecos_set), ";".join(item[1]), # for all
                            "", "", ""# Only for large scale study data
                            ])


    # find PTMS in EBI API            
    url_ebi = f"https://www.ebi.ac.uk/proteins/api/proteomics/ptm/{accession}"

            
    try:
        response_ebi = requests.get(url_ebi, headers=headers)
        response_ebi.raise_for_status()
    except requests.exceptions.HTTPError as e:
        # logging.info(f"Entry with accession {accession} not found on EBI API.")
        end_time = time.time()
        elapsed = end_time - start_time

        logging.info(f"Not on EBI API. Time elapsed for {accession}: {elapsed:.2f} seconds.")
        return df_rows, seq_rows

    data_ebi = response_ebi.json()

    # to deduplicate the PTMs from EBI
    ptms_ebi = {}

    for feature in data_ebi.get("features", []):
        feature_type = feature.get("type", "")
        if feature_type not in ["PROTEOMICS_PTM"]:
            continue

        # peptide start position
        pep_start = int(feature.get("begin", 0))

        #eco codes
        ecos = set()
        for evidence in feature.get("evidences", []):
            eco = evidence.get("code","")
            if eco not in valid_ECO_codes:
                continue
            ecos.add(eco)

        # PTM description
        for ptm in feature.get("ptms", []):
                    
            #description
            desc = ptm.get("name", "").replace(" ", "__")

            # get modified amino acid position 
            pos_in_peptide = int(ptm.get("position", 0))
            position = pep_start + pos_in_peptide - 1

            lss_databases = set(ptm.get("sources", []))

            #ptmxchange ids and sources and confidence scores
            x_ids = set()
            pubmed_ids = set()
            conf_scores = set()
            for reference in ptm.get("dbReferences", []):
                x_ids.add(reference.get("id", ""))
                pubmed_ids.add(reference.get("properties", {}).get("Pubmed ID", ""))
                conf_scores.add(reference.get("properties", {}).get("Confidence score", ""))
                    
            # deduplicates based on position + description
            ptm_id = f"{position}|{desc}"
            ptm_desc = [ecos, lss_databases, x_ids, pubmed_ids, conf_scores]

            if ptm_id not in ptms_ebi:
                ptms_ebi[ptm_id] = ptm_desc
            else:
                ptms_ebi[ptm_id][0].update(ecos)
                ptms_ebi[ptm_id][1].update(lss_databases)
                ptms_ebi[ptm_id][2].update(x_ids)
                ptms_ebi[ptm_id][3].update(pubmed_ids)
                ptms_ebi[ptm_id][4].update(conf_scores)

    # appends to dataframe
    for key, item in ptms_ebi.items():
        ecos_set = item[0]
        if not ecos_set:
            continue
        if item[3] == {} or item[3] == {""}:
            pubmed_ids_str = ""
        else:
            pubmed_ids_str = f"PubMed:{';'.join(item[3])}"
        df_rows.append([accession, key.split("|")[0], length, entry_type,
                            protein_name, "PROTEOMICS__PTM", key.split("|")[1], organism_name,
                            ";".join(ecos_set), pubmed_ids_str, # for all
                            ";".join(item[1]), ";".join(item[2]), ":".join(item[4]) # Only for large scale study data
                            ]) 
                
                        
    end_time = time.time()
    elapsed = end_time - start_time

    logging.info(f"Time elapsed for {accession}: {elapsed:.2f} seconds.")

    return df_rows, seq_rows

def get_uniprot_records(data_dir, workers=4):
    """
    Fetch records using threading."""

    # temporary file output path
    out_info = f"{data_dir}/deduplicated_records.csv.tmp"
    out_seq = f"{data_dir}/sequences.csv.tmp"

    # REST URL
    base_url_rest = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url_rest = f"{base_url_rest}?query=existence:1&format=json&size=500"

    # for pagination purposes
    tally = 0

    # to store data
    dataframe = []
    # to store sequences
    sequences = []

    logging.info(f"Fetching entries from UniProt, please wait.")

    while url_rest:
        response = requests.get(url_rest, headers=headers)
        response.raise_for_status()
        data = response.json()

        results = data.get("results", [])

        with ThreadPoolExecutor(max_workers=workers) as executor:  
            futures = [
                executor.submit(process_entry, entry, headers)
                for entry in results
            ]
            for fut in as_completed(futures):
                df_rows, seq_rows = fut.result()
                dataframe.extend(df_rows)
                sequences.extend(seq_rows)

        # writes to temporary csv to avoid data loss
        df = pd.DataFrame(dataframe, columns=columns)
        seq_df = pd.DataFrame(sequences, columns=seq_columns)

        if not os.path.exists(out_info) or os.path.getsize(out_info) == 0:
            df.to_csv(out_info, index=False, mode="w")        # includes header
        else:
            df.to_csv(out_info, index=False, mode="a", header=False)

        # Free memory for next page
        dataframe = []
        del df

        if not os.path.exists(out_seq) or os.path.getsize(out_seq) == 0:
            seq_df.to_csv(out_seq, index=False, mode="w")        # includes header
        else:
            seq_df.to_csv(out_seq, index=False, mode="a", header=False)
        # Free memory for next page
        sequences = []
        del seq_df

        # pagination
        tally += len(data["results"])
        total = response.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)

        url_rest = response.links.get("next", {}).get("url")
        if tally >= int(total):
            break

    logging.info(f"Entry download is done. {os.path.getsize(out_info)} records stored at {out_info}.") 
    logging.info(f"{os.path.getsize(out_seq)} sequences stored at {out_seq}.") 

if __name__ == "__main__":   
    import argparse          

    parser = argparse.ArgumentParser(description="Fetch UniProt PTM records.")
    parser.add_argument("data_dir", help="Output directory")  
    parser.add_argument("--threads", type=int, default=8,
                        help="Number of threads for parallel requests")  

    args = parser.parse_args()

    get_uniprot_records(args.data_dir, workers=args.threads) 
    