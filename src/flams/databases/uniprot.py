#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ieva-balta, majocava, naaattella
"""

import sys
import os
import json
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import requests
import logging
import re

from concurrent.futures import ThreadPoolExecutor, as_completed

from . import setup

""" uniprot
This script downloads the different contents of the UniProt database, and transforms them into a fasta format.
Script developed to work with UniProt database release 2025_04.
"""

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

#valid eco codes
valid_ECO_codes = [
    "ECO:0000269", #(experimental evidence used in manual assertion)​
    "ECO:0000314", #(direct assay evidence used in manual assertion)​
    "ECO:0007744", #(combinatorial computational and experimental evidence used in manual assertion)​
    "ECO:0007829", #(combinatorial computational and experimental evidence used in automatic assertion)​
]

# columns for the records files
info_columns = ["Accession", "Position", 
                            "Feature_type", "Description", 
                            "ECO_codes", "Sources", # for all
                            "LSS_Database", "LSS_IDs", "LSS_Confidence_scores"# Only for large scale study data (from EBI API)
                            ]

# dtypes for the records files
info_dtypes = {
    "Accession": "string",
    "Position": "Int64",
    "Feature_type": "string",
    "Description": "string",
    "ECO_codes": "string",
    "Sources": "string",
    "LSS_Database": "string",
    "LSS_IDs": "string",
    "LSS_Confidence_scores": "string",
}

# columns for the general sequence  info file
general_columns = ["Accession", "Length", "Entry_type", "Protein", "Organism","Sequence"]

# dtypes for the general sequence info file
general_dtypes = {
    "Accession": "string",
    "Length": "Int64",
    "Entry_type": "string",
    "Protein": "string",
    "Organism": "string",
    "Sequence": "string",
}

def get_fasta(PTM_modification_dict, data_dir):
    """

    This function downloads all entries of the UniProt database, using two different APIs (REST and EBI).
    It sorts them based on regular expressions and saves them in fasta format in $data_dir.
    Additionally, the function logs unclassified entries. 

    Parameters
    ----------
    PTM_modification_dict: Dict[]
        Dictionary of PTM types and info associated with it, such as, version, RegEx for sorting and list of amino acids.
    data_dir: str
        Location where output files will be stored
    """

    # temporary file that stores merged version of rest and ebi outputs
    # the file could be removed later
    records_file = os.path.join(data_dir, f"_deduplicated_records-{setup.version}.csv.tmp")

    # checks if records file already exists
    if not os.path.exists(records_file) or os.path.getsize(records_file) == 0:
        logging.info(f"Records file not found. Proceeding to download entries from UniProt.")
        # downloads entries from rest api and stores accession codes for ebi api
        accessions = download_with_rest_api(data_dir, version=setup.version)
        # downloads entries from ebi api using accession codes from rest api
        download_with_ebi_api(data_dir, accessions, version=setup.version)
        # merges both files into one
        merge_rest_ebi_outputs(data_dir)
        df = pd.read_csv(records_file, header=0, low_memory=False)
        df = df.astype(info_dtypes)
    else:
        df = pd.read_csv(records_file, header=0, low_memory=False)
        df = df.astype(info_dtypes)

    # classifies records into PTM types and stores fasta files per modification type
    logging.info(f"Sorting records into modification types and creating fasta files.")
    classified = sort_uniprot_records(df, PTM_modification_dict, data_dir)

    # finds unclassified records 
    unclassified = pd.concat([df, classified]).drop_duplicates(keep=False)

    ### for development purposes - can remove later
    # stores a list of classified and unclassified records
    classified.to_csv(f"{data_dir}/classified-{setup.version}.csv", index=False)

    # stores unclassified records
    # this is what logs unclassified entries basically
    unclassified.to_csv(f"{data_dir}/unclassified-{setup.version}.csv", index=False)
    logging.info(f"CSV file for unclassified entries was created and stored at {data_dir}/unclassified-{setup.version}.csv.")


def sort_uniprot_records(uniprot_records, PTM_modification_dict, data_dir):
    """

    This function groups uniprot records by PTM type, saves the entries in fasta format.
    Additionally, it returns all classified records as a dataframe.

    Parameters
    ----------
    uniprot_records: pd.Dataframe
        Pandas dataframe of uniprot records
    PTM_modification_dict: Dict[]
        Dictionary of PTM types and info associated with it, such as, version, RegEx for sorting and list of amino acids.
    data_dir: str
        Location where output files will be stored

    """
    # make new dataframe
    classified_records = pd.DataFrame({
        col: pd.Series(dtype=info_dtypes[col]) 
        for col in info_columnsgit 
    })

    # parses through modification dictionary
    for modification, m_type in PTM_modification_dict.items():
       
        # finds all records matching the specific RegEx
        df_mod_type = find_modification_type(uniprot_records, m_type, data_dir)

        # skips if no records found for this modification type
        if df_mod_type.empty:
            continue
        
        # adds the records of this modification to the dataframe of all classified records
        classified_records = pd.concat([classified_records, df_mod_type], ignore_index=True)
    
    return classified_records
    

def find_modification_type(uniprot_records, PTM_modification_type, data_dir):
    """

    This function finds and returns all records where the PTM descriptions match the specified regular expressions.

    Parameters
    ----------
    uniprot_records: pd.Dataframe
        Pandas dataframe of uniprot records
    PTM_modification_type: ModificationType
        ModificationType object containing info about the modification type
    data_dir: str
        Location where output files will be stored

    """
    # retrieves modification type info
    m_type = PTM_modification_type.type
    m_version = PTM_modification_type.version
    m_db = PTM_modification_type.dbs[0]
    regex = m_db.descriptor

    # combines all regex patterns into one (using OR)
    combined_pattern = "|".join(f"(?:{p.lower()})" for p in regex)

    # matches RegEx's
    description = uniprot_records["Description"].astype(str).str.lower()
    mask = description.str.contains(combined_pattern, na=False, regex=True, flags=re.IGNORECASE)

    # filters the input dataframe
    df_filtered = uniprot_records.loc[mask].copy()

    # skips empty matches - only useful when downloading a new version of uniprot or changing the modification types in setup.MODIFICATIONS
    if df_filtered.empty:
        logging.info(f"No entries for modification type {m_type} were found.")
        return df_filtered
    
    # passes data to fasta conversion
    df_to_fasta(df_filtered, m_type, m_version, data_dir)

    ### for development purposes - can remove later
    # stores records per modification type
    df_filtered.to_csv(f"{data_dir}/{m_type}-{m_version}.csv", index=False)

    return df_filtered

def df_to_fasta(PTM_type_df, m_type, m_version, data_dir):
    """

    This function writes the PTM records into a fasta format.

    Parameters
    ----------
    PTM_type_df: pd.Dataframe
        Pandas dataframe of uniprot records per modification type
    m_type: str
        Modification type name
    m_version: str
        Database version
    data_dir: str
        Location where output files will be stored
    """
    # reads sequence file containing general information on each entry
    sequence_file = os.path.join(data_dir, f"_sequences-{setup.version}.csv.tmp")
    seq_df = pd.read_csv(sequence_file, header=0, low_memory=False)
    seq_df = seq_df.astype(general_dtypes)

    # adds sequences and general information to PTM dataframe
    PTM_df = PTM_type_df.merge(seq_df, how="left", on="Accession")

    # counts number of records written
    tally = 0

    # writes to a fasta file
    with open(os.path.join(data_dir, f"{m_type}-{m_version}.fasta"), "a") as f:

        for _, row in PTM_df.iterrows():

            # replaces the NoneType  with N/A (needed for run_blast.py)
            protein = row["Protein"]
            if pd.isna(protein) or protein == "":
                protein = "N/A"

            sources = row["Sources"]
            if pd.isna(sources) or sources == "":
                sources = "N/A"

            lss_database = row["LSS_Database"]
            if pd.isna(lss_database) or lss_database == "":
                lss_database = "N/A"

            lss_ids = row["LSS_IDs"]
            if pd.isna(lss_ids) or lss_ids == "":
                lss_ids = "N/A"

            lss_confidence_scores = row["LSS_Confidence_scores"]
            if pd.isna(lss_confidence_scores) or lss_confidence_scores == "":
                lss_confidence_scores = "N/A"
            
            # adds unique id instead of accession (needed for making BLAST databases)
            unique_id = f"{row['Accession']}_{tally}"
            tally += 1
            
            # formats the sequence
            seq = Seq(row["Sequence"])

            # sequence id
            id = f'{unique_id}|{str(row["Position"])}|{row["Length"]}|{row["Entry_type"]}'
            # formats the record
            record = SeqRecord(
                    seq,
                    id=id,
                    description=f'{protein}|{row["Description"]}|{row["Organism"]} [{row["ECO_codes"]}|{sources}|{lss_database}|{lss_ids}|{lss_confidence_scores}]',
                )
            # writes to fasta
            SeqIO.write(record, f, "fasta")
        logging.info(f"Fasta file with {tally} records for modifictation '{m_type}' created and stored at '{data_dir}'.")

def merge_rest_ebi_outputs(data_dir):
    """
    
    This function merges the temporary files from REST and EBI APIs into one deduplicated records file.

    Parameters
    ----------
    data_dir: str
        Location where output files will be stored
    """
    # file paths
    rest_file = os.path.join(data_dir, f"_deduplicated_rest_api_records-{setup.version}.csv.tmp")
    ebi_file = os.path.join(data_dir, f"_deduplicated_ebi_api_records-{setup.version}.csv.tmp")
    merge_file = os.path.join(data_dir, f"_deduplicated_records-{setup.version}.csv.tmp")

    # checks if both files exist
    if not os.path.exists(rest_file) or os.path.getsize(rest_file) == 0:
        logging.info(f"No REST API records found to merge.")
        return
    if not os.path.exists(ebi_file) or os.path.getsize(ebi_file) == 0:
        logging.info(f"No EBI API records found to merge.")
        return

    # reads rest and ebi files
    df_rest = pd.read_csv(rest_file, header=0, low_memory=False)
    df_rest = df_rest.astype(info_dtypes)
    df_ebi = pd.read_csv(ebi_file, header=0, low_memory=False)
    df_ebi = df_ebi.astype(info_dtypes)

    # merges rest and ebi dataframes
    merged_df = pd.concat([df_rest, df_ebi], ignore_index=True)

    # writes merged dataframe to file
    merged_df.to_csv(merge_file, index=False)
    logging.info(f"Merged records file created and stored at {merge_file}.")

def download_with_ebi_api(data_dir, accession_list, version=None, threads=40):
    """
    
    This funnction downloads PTM entries from the EBI API using a list of accession codes from the REST API.
    The EBI API only fetches PTM entries from large scale studies, complementing the data from the REST API.

    Parameters
    ----------
    data_dir: str
        Location where output files will be stored
    accession_list: list
        List of accession codes to download
    version: str, optional
        Version identifier for the download (default is None)
    threads: int, optional
        Number of threads to use for downloading (default is 40)
        Used to speed up the download process by parallelizing requests.
    """
    # if version not specified, use setup version
    if version is None:
        version = setup.version

    # temporary file output paths
    out_info = f"{data_dir}/_deduplicated_ebi_api_records-{version}.csv.tmp"
    processed_path = f"{data_dir}/_processed_ebi_accessions-{version}.txt"

    # EBI API URL
    base_url_ebi = "https://www.ebi.ac.uk/proteins/api/proteomics/ptm"
    headers = {"Accept": "application/json"}

    logging.info(f"Starting entry download using EBI API. Total accessions to go through: {len(accession_list)}")

    # stores processed accessions
    processed = set()
    if os.path.isfile(processed_path):
        with open(processed_path) as f:
            processed = {line.strip() for line in f}

    # only process unprocessed accessions
    to_process = [acc for acc in accession_list if acc not in processed]
    logging.info(f"Total remaining accessions: {len(to_process)}")

    # create batches of accessions
    batch_size = 100
    batches = list(batch_accessions(to_process, batch_size))
    logging.info(f"Total batches to process: {len(batches)}")

    # uses ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=threads) as ex:

        # submits each batch to the thread pool
        futures = {ex.submit(fetch_ptm_for_batch, batch, base_url_ebi, headers): batch for batch in batches}

        # counts processed accessions
        new_batch_processed = 0

        for f in as_completed(futures):
            current_batch = futures[f]
            new_batch_processed += len(current_batch)
            try:
                rows, processed_accs = f.result()
                
                if rows:
                    write_ebi_output(rows, out_info)
                    
                logging.info(
                    f"Batch processed ({new_batch_processed} accs): "
                    f"{len(rows)} PTM sites found."
                )

                # updates processed.txt with all accessions in the batch
                with open(processed_path, "a") as p:
                    for acc in current_batch:
                        p.write(acc + "\n")

            except Exception as e:
                # logs an error for the batch but do NOT mark as processed
                logging.error(
                    f"ERROR processing batch of {len(current_batch)} accessions "
                    f"starting with {current_batch[0]}: {e}"
                )

    logging.info("PTM download using EBI API is done.")

def write_ebi_output(all_df_rows, out_info):
    """
    
    This function writes the EBI API PTM data to a CSV file.

    Parameters
    ----------
    all_df_rows: list
        List of rows to write to the CSV file
    out_info: str
        Output file path for the CSV file
    """
    df = pd.DataFrame(all_df_rows, columns=info_columns)
    df = df.astype(info_dtypes)
    if not os.path.exists(out_info) or os.path.getsize(out_info) == 0:
        df.to_csv(out_info, index=False, mode="w")        # includes header
    else:
        df.to_csv(out_info, index=False, mode="a", header=False)
def fetch_ptm_for_batch(batch_accessions, base_url_ebi, headers):
    """
    
    This function fetches PTM data from the EBI API for a batch of accession codes.

    Parameters
    ----------
    batch_accessions: list
        List of accession codes
    base_url_ebi: str
        Base URL for the EBI API
    headers: dict
        Headers for the API request
    """
    # constructs the URL for the batch request
    acc_list = ",".join(batch_accessions)
    url_ebi = f"{base_url_ebi}?accession={acc_list}&size=100"

    # sends the GET request to the EBI API
    response = requests.get(url_ebi, headers=headers)
    response.raise_for_status()
    data = response.json()

    # initializes lists to store results
    all_df_rows = []
    processed_accessions = set()

    # parses through each entry in batch
    for entry in data:
        rows = process_ebi_entry(entry)
        all_df_rows.extend(rows)

        if "accession" in entry:
            processed_accessions.add(entry["accession"])

    return all_df_rows, list(processed_accessions)


def process_ebi_entry(entry):
    """
    
    This function processes a single entry from the EBI API and extracts PTM data.

    Parameters
    ----------
    entry: dict
        A single entry from the EBI API response
    """

    # to store dataframe rows
    df_rows = []

    # to deduplicate the PTM information from EBI
    ptms_ebi = {}

    accession = entry.get("accession", "")

    # parses through features
    for feature in entry.get("features", []):
        feature_type = feature.get("type", "")
        if feature_type not in ["PROTEOMICS_PTM"]:
            continue

        # peptide start position
        pep_start = int(feature.get("begin", 0))

        # gets ECO codes
        ecos = set()
        for evidence in feature.get("evidences", []):
            eco = evidence.get("code","")
            # skips if ECO codes are not valid
            if eco not in valid_ECO_codes:
                continue
            ecos.add(eco)

        # PTM description
        for ptm in feature.get("ptms", []):
                    
            # description
            desc = ptm.get("name", "").replace(" ", "__")

            # gets modified amino acid position 
            pos_in_peptide = int(ptm.get("position", 0))
            position = pep_start + pos_in_peptide - 1

            # large scale study databases - PTMXchange or PRIDE
            lss_databases = set(ptm.get("sources", []))

            # LSS ids (x_ids), sources and confidence scores
            x_ids = set()
            pubmed_ids = set()
            conf_scores = set()
            for reference in ptm.get("dbReferences", []):
                x_ids.add(reference.get("id", "").replace(" ", "__"))
                pubmed_ids.add(reference.get("properties", {}).get("Pubmed ID", "").replace(" ", "__"))
                conf_scores.add(reference.get("properties", {}).get("Confidence score", "").replace(" ", "__"))
                    
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
        # skips if no valid ECO codes
        if not ecos_set:
            continue
        # formats pubmed ids if there aren't any
        if item[3] == {} or item[3] == {""}:
            pubmed_ids_str = ""
        else:
            pubmed_ids_str = f"PubMed:{';'.join(item[3])}"
        df_rows.append([accession, key.split("|")[0],
                            "PROTEOMICS__PTM", key.split("|")[1],
                            ";".join(ecos_set), pubmed_ids_str, # for all
                            ";".join(item[1]), ";".join(item[2]), ":".join(item[4]) # Only for large scale study data
                            ]) 

    return df_rows

def batch_accessions(accession_list, batch_size):
    """

    Batches a list of accession numbers into smaller lists of a specified size.
    The last batch might not be the exact size.

    Parameters
    ----------
    accession_list: list
        List of accession numbers to be batched
    batch_size: int
        Size of each batch
    """
    for i in range(0, len(accession_list), batch_size):
        yield accession_list[i:i + batch_size]

def download_with_rest_api(data_dir, version=None):
    """
    
    This function downloads PTM entries from the UniProt REST API.
    The REST API fetches general PTM entries, excluding large scale study data.

    Parameters
    ----------
    data_dir: str
        Location where output files will be stored
    version: str, optional
        Version identifier for the download (default is None)
    """
    # if version not specified, use setup version
    if version is None:
        version = setup.version
        
    # temporary file output path
    out_info = f"{data_dir}/_deduplicated_rest_api_records-{version}.csv.tmp"
    out_seq = f"{data_dir}/_sequences-{version}.csv.tmp"

    # REST URL
    base_url_rest = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url_rest = f"{base_url_rest}?query=existence:1&includeIsoform=true&format=json&size=500"

    # for pagination purposes
    tally = 0

    # to store data
    dataframe = []
    # to store general info separately (useful for EBI API)
    general_info = []

    # to store accession code for EBI API
    accession_list = []

    logging.info(f"Fetching entries from UniProt, please wait.")

    while url_rest:
        response = requests.get(url_rest, headers=headers)
        response.raise_for_status()
        data = response.json()

        results = data.get("results", [])

        # parse through results
        for entry in results:
            df_rows = []
            gen_rows = []

            # gets general information
            accession = entry.get("primaryAccession")
            sequence = entry.get("sequence", {}).get("value", "")
            length = len(sequence)
            name = entry.get("genes", [{}])[0].get("geneName", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")

            # database
            entry_type = entry.get("entryType", "")
            if "Swiss-Prot" in entry_type:
                entry_type = "Swiss-Prot"
            if "TrEMBL" in entry_type:
                entry_type = "TrEMBL"
            
            # removes white space
            entry_type = entry_type.replace(" ", "__")
            protein_name = f"{name}".replace(" ","__")
            organism_name = f"{organism}".replace(" ","__")

            # stores accession for EBI API
            accession_list.append(accession)
            # stores sequence separately
            gen_rows.append([accession, length, entry_type, protein_name, organism_name, sequence])
                
            # finds PTM sites

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

                # removes the notes from the description
                if ";" in desc:
                    desc = desc.split(";")[0].strip()

                # special case - "removed"
                if desc.lower() == "removed":
                    continue

                # special case - microbial infection
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

                # in case of different start and end positions adds both positions
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
                # skips if no valid ECO codes
                if not ecos_set:
                    continue
                df_rows.append([accession, key.split("|")[0], 
                            item[2], key.split("|")[1], 
                            ";".join(ecos_set), ";".join(item[1]), # for all
                            "", "", "" # Only for large scale study data
                            ])

            dataframe.extend(df_rows)
            general_info.extend(gen_rows)

        # writes to a temporary csv to avoid data loss
        df = pd.DataFrame(dataframe, columns=info_columns)
        df = df.astype(info_dtypes)
        general_df = pd.DataFrame(general_info, columns=general_columns)
        general_df = general_df.astype(general_dtypes)

        if not os.path.exists(out_info) or os.path.getsize(out_info) == 0:
            df.to_csv(out_info, index=False, mode="w")        # includes header
        else:
            df.to_csv(out_info, index=False, mode="a", header=False)

        # free memory for next page
        dataframe = []
        del df

        if not os.path.exists(out_seq) or os.path.getsize(out_seq) == 0:
            general_df.to_csv(out_seq, index=False, mode="w")        # includes header
        else:
            general_df.to_csv(out_seq, index=False, mode="a", header=False)
            
        # free memory for next page
        general_info = []
        del general_df

        # pagination
        tally += len(data["results"])
        total = response.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results using REST API", tally, total)

        url_rest = response.links.get("next", {}).get("url")
        if tally >= int(total):
            break

    logging.info(f"Entry download using REST API is done. {os.path.getsize(out_info)} records stored at {out_info}.") 
    logging.info(f"General information for entries stored at {out_seq}.")
    
    # returns accession list for EBI API
    return accession_list

