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

### for development purposes
#  to avoid running download all the time
# df = pd.read_csv("/data/leuven/368/vsc36826/flams/flams2/dbs/deduplicated_records-2.0.csv", header=0)

def get_fasta(PTM_modification_dict, data_dir):
    """
    make fasta files
    log unclassified records

    This function downloads all entries of the UniProt database,
    and sorts them based on regular expressions and saves them in fasta format in $data_dir.
    Additionally, the function logs unclassified entries. 

    Parameters
    ----------
    PTM_modification_dict: Dict[]
        Dictionary of PTM types and info associated with it, such as, version, RegEx for sorting and list of amino acids.
    data_dir: str
        Location where output files will be stored
    """

    ### for development purposes
    # all_records = df
    # downloads and deduplicates all records
    # all_records = get_uniprot_records(data_dir)

    ### checks for the file, 
    ### think about how to integrate the downlaod into the pipeline properly
    records_file = os.path.join(data_dir, f"_deduplicated_records-{setup.version}.csv.tmp")

    if not os.path.exists(records_file) or os.path.getsize(records_file) == 0:
        logging.error(f"Records file not found. Please make sure to download UniProt records first, using 'download_uniprot_records.py' script. Exiting FLAMS...")
        sys.exit()
    else:
        df = pd.read_csv(records_file, header=0)

    sequence_file = os.path.join(data_dir, f"_sequences-{setup.version}.csv.tmp")

    if not os.path.exists(sequence_file) or os.path.getsize(sequence_file) == 0:
        logging.error(f"Sequence file not found. Please make sure to download UniProt records first, using 'download_uniprot_records.py' script. Exiting FLAMS...")
        sys.exit()


    # classifies records into PTM types and stores fasta files per modification type
    classified = sort_uniprot_records(df, PTM_modification_dict, data_dir)

    # finds unclassified records 
    unclassified = pd.concat([df, classified]).drop_duplicates(keep=False)

    ### for development purposes - can remove later
    # stores a list of classified and unclassified records
    classified.to_csv(f"{data_dir}/classified-{setup.version}.csv", index=False)

    #stores unclassified records
    unclassified.to_csv(f"{data_dir}/unclassified-{setup.version}.csv", index=False)

    logging.info(f"CSV file for unclassified entries was created and stored at {data_dir}/unclassified-{setup.version}.csv.")
    # # outputs unclassified records in a fasta format
    # fasta_records_unclassified = df_to_fasta(unclassified)
    # with open(f"{data_dir}/unclassified-{setup.version}.fasta", "w", encoding="UTF-8") as out:
    #     SeqIO.write(fasta_records_unclassified, out, "fasta")

    # logging.info(f"Fasta file for {len(unclassified)} unclassified entries was created and stored at {data_dir}/unclassified-{setup.version}.fasta.")


def sort_uniprot_records(uniprot_records, PTM_modification_dict, data_dir):
    """

    This function groups uniprot records by PTM type, saves the entries in fasta format.
    Additionally it returns all classified records as a dataframe.


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
    classified_records = pd.DataFrame(columns = ["Accession", "Position", "Length", "Entry_type",
                                                    "Protein", "Feature_type", "Description", "Organism",
                                                    "ECO_codes", "Sources", 
                                                    "LSS_Database", "LSS_IDs", "LSS_Confidence_scores"])

    # parses through modification dictionary
    for modification, m_type in PTM_modification_dict.items():
        # # gets the RegEx descriptors for the modification type
        # m_db = m_type.dbs[0]
        # regex = m_db.descriptor
        # finds all records matching the specific RegEx
        df_mod_type = find_modification_type(uniprot_records, m_type, data_dir)

        if df_mod_type.empty:
            continue

        # ### for development purposes - can remove later
        # # stores records per modification type
        # df_mod_type.to_csv(f"{data_dir}/{modification}-{m_type.version}.csv", index=False)
        
        # adds the records of this modification to the dataframe of all classified records
        classified_records = pd.concat([classified_records, df_mod_type], ignore_index=True)

        # # outputs a fasta file per modification type
        # fasta_records = df_to_fasta(df_mod_type)
        # with open(f"{data_dir}/{modification}-{m_type.version}.fasta", "w", encoding="UTF-8") as out:
        #     SeqIO.write(fasta_records, out, "fasta")
        
        # logging.info(f"Fasta file with {len(df_mod_type)} records for modifictation '{modification}' created and stored at {data_dir}.")
    
    return classified_records
    

def find_modification_type(uniprot_records, PTM_modification_type, data_dir):
    """

    This function finds and returns all records where the PTM descriptions matches the specified regular expressions.

    Parameters
    ----------
    uniprot_records: pd.Dataframe
        Pandas dataframe of uniprot records
    PTM_regex_list: List[str]
        List of regular expressions that should be matched

    """
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

    #skip empty mathces - only useful when downloading a new version of uniprot or changing the modification type sin setup.MODIFICATIONS
    if df_filtered.empty:
        logging.info(f"No entries for modification type {m_type} were found.")
        return df_filtered
    
    # pass to fasta conversion
    df_to_fasta(df_filtered, m_type, m_version, data_dir)

    # for development purposes
    # stores records per modification type
    df_filtered.to_csv(f"{data_dir}/{m_type}-{m_version}.csv", index=False)

    return df_filtered

def df_to_fasta(PTM_type_df, m_type, m_version, data_dir):
    """

    This function formats a pandas dataframe into a fasta records list.

    Parameters
    ----------
    PTM_type_df: pd.Dataframe
        Pandas dataframe of uniprot records
    """
    sequence_file = os.path.join(data_dir, f"_sequences-{setup.version}.csv.tmp")
    seq_df = pd.read_csv(sequence_file, header=0)

    # adds sequences to PTM dataframe
    PTM_df = PTM_type_df.merge(seq_df, how="left", on="Accession")

    suffix = 0
    # change name
    with open(os.path.join(data_dir, f"{m_type}-{m_version}.fasta"), "a") as f:
        for _, row in PTM_df.iterrows():
            # replaces the NoneType  with N/A (needed for run_blast.py)
            protein = row["Protein"] or "N/A"
            sources = row["Sources"] or "N/A"
            lss_database = row["LSS_Database"] or "N/A"
            lss_ids = row["LSS_IDs"] or "N/A"
            lss_confidence_scores = row["LSS_Confidence_scores"] or "N/A"
            
            # adds unique id instead of accession (needed for making BLAST databases)
            unique_id = f"{row['Accession']}_{suffix}"
            suffix += 1
            
            #sequence
            seq = Seq(row["Sequence"])

            # sequence id
            id = f'{unique_id}|{str(row["Position"])}|{row["Length"]}|{row["Entry_type"]}'
            record = SeqRecord(
                    seq,
                    id=id,
                    description=f'{protein}|{row["Description"]}|{row["Organism"]} [{row["ECO_codes"]}|{sources}|{lss_database}|{lss_ids}|{lss_confidence_scores}]',
                )

            SeqIO.write(record, f, "fasta")
        logging.info(f"Fasta file with {suffix} records for modifictation '{m_type}' created and stored at {data_dir}.")    

    # # parses through the dataframe
    # for idx, row in PTM_type_df.iterrows():
    #     # replaces the NoneType protein names with N/A (needed for run_blast.py)
    #     protein = row["Protein"]
    #     if not protein:
    #         protein = "N/A"
        
    #     seq = Seq(row["Sequence"])
    #     # adds unique id instead of accession (needed for making BLAST databases)
    #     id = f'{row["Unique_ID"]}|{str(row["Position"])}|{row["Length"]}|{row["Entry_type"]}'
    #     rec = SeqRecord(
    #                 seq,
    #                 id=id,
    #                 description=f'{protein}|{row["Description"]}|{row["Organism"]} [{row["ECO_codes"]}|{row["Sources"]}|{row["Source_ids"]}]',
    #             )

    #     fasta_records.append(rec)

    # return fasta_records

# def deduplicate_records(uniprot_records):
#     """

#     This function takes a pandas dataframe of uniprot records and deduplicates it
#     based on the combination of Accession, PTM position, and PTM description.
#     Additionally, it concatinates the ECO codes and source information from the duplicated entries.

#     Parameters
#     ----------
#     uniprot_records: pd.Dataframes
#         Pandas dataframe of uniprot records

#     """

#     # columns to group by
#     group_cols = ["Accession", "Position", "Description"]

#     # aggregated dataframes
#     collapsed = (
#         uniprot_records.groupby(group_cols, as_index=False)
#         .agg({
#             "Unique_ID": "first",
#             "Accession": "first",
#             "Position": "first",
#             "Length": "first",
#             "Entry_type": "first",
#             "Protein": "first",
#             "Feature_type": "first",
#             "Description": "first",
#             "Organism": "first",
#             # combine these fields by joining unique non-null values
#             "ECO_codes": lambda x: ";".join(sorted(set(filter(None, x)))),
#             "Sources": lambda x: ";".join(sorted(set(filter(None, x)))),
#             "Source_ids": lambda x: ";".join(sorted(set(filter(None, x)))),
#             "Sequence": "first"
#         })
#     )

#     return collapsed
    
# def get_uniprot_records(data_dir):
#     """
#     fetches all uniprot ptm records
#     outputs a dataframe of all records + info

#     This function fetches all relevant uniprot entries (existence:1 and valid_eco_codes) and stores them in a pandas dataframe.

#     Parameters
#     ----------
#     data_dir: str
#         Location for storing data

#     """
#     logging.info(f"Fetching entries from UniProt, please wait.")

#     base_url = "https://rest.uniprot.org/uniprotkb/search"
#     headers = {"Accept": "application/json"}
#     url = f"{base_url}?query=existence:1&format=json&size=500"

#     # for pagination purposes
#     tally = 0

#     # to store data
#     dataframe = []

#     # to keep track of modified residues all together
#     tally_per_feature = 0

#     while url:
#         response = requests.get(url, headers=headers)
#         response.raise_for_status()
#         data = response.json()

#        # parse through results
#         for entry in data.get("results", []):
#             accession = entry.get("primaryAccession")
#             sequence = entry.get("sequence", {}).get("value", "")
#             name = entry.get("genes", [{}])[0].get("geneName", {}).get("value", "")
#             organism = entry.get("organism", {}).get("scientificName", "")

#             # storing database
#             entry_type = entry.get("entryType", "")
#             if "Swiss-Prot" in entry_type:
#                 entry_type = "Swiss-Prot"
#             if "TrEMBL" in entry_type:
#                 entry_type = "TrEMBL"

#             # removes white space
#             entry_type = entry_type.replace(" ", "__")
#             protein_name = f"{name}".replace(" ","__")
#             organism_name = f"{organism}".replace(" ","__")
#             length = len(sequence)
            
#             # finds PTM sites from REST
#             # finds relevant feature types
#             for feature in entry.get("features", []):
#                 feature_type = feature.get("type")
#                 if feature.get("type") not in ["Modified residue",
#                             # "Chain",
#                             "Modified residue (large scale data)", 
#                             "Lipidation", 
#                             "Glycosylation",
#                             "Disulfide bond",
#                             "Cross-link"]:
#                     continue

#                 # PTM description
#                 desc = feature.get("description", "")

#                 #special case  - disulfide bonds
#                 if feature_type == "Disulfide bond":
#                     desc = "Disulfide bond"

#                 # remove the notes from the description
#                 if ";" in desc:
#                     desc = desc.split(";")[0].strip()

#                 #special case - "removed"
#                 if desc.lower() == "removed":
#                     continue

#                 #special case - microbial infection
#                 if "microbial infection" in desc.lower():
#                     desc = desc.replace("(Microbial infection) ", "")

#                 # gets positions for start and end (matters for crosslinks, bond formation, to record both amino acids)
#                 pos_start = feature.get("location", {}).get("start", {}).get("value", "N/A")
#                 pos_end = feature.get("location", {}).get("end", {}).get("value", "N/A")

#                 # gets evidence ECO|source|ids
#                 ECOs = []
#                 sources = []
#                 ids = []
#                 for ev in feature.get("evidences", []):
#                     eco = ev.get("evidenceCode", "")
#                     # skips if the evidnece is from the valid eco codes
#                     if eco not in valid_ECO_codes:
#                         continue
#                     source = ev.get("source", "")
#                     source_id = str(ev.get("id", ""))
#                     ECOs.append(eco)
#                     sources.append(source)
#                     ids.append(source_id)
#                 ECO_str = ";".join(ECOs) 
#                 # skips if there are no ECOs
#                 if not ECOs:
#                     continue 
#                 source_str = ";".join(sources) if sources else ""
#                 ids_str = ";".join(ids) if ids else ""

                

#                 # removes white space
#                 feature_type_name = f"{feature_type}".replace(" ", "__")
#                 desc = f"{desc}".replace(" ","__")
                
#                 # adds ID to differentiate
#                 unique_id = f"{accession}_{tally_per_feature}"
#                 tally_per_feature += 1

#                 # appends to dataframe
#                 # add with start position
#                 dataframe.append([unique_id, accession, pos_start, length, entry_type, 
#                                     protein_name, feature_type_name, desc, organism_name, 
#                                     ECO_str, source_str, ids_str, sequence])

#                 # adds ID to differentiate
#                 unique_id = f"{accession}_{tally_per_feature}"
#                 tally_per_feature += 1
                
#                 # adds the end position
#                 dataframe.append([unique_id, accession, pos_end, length, entry_type, 
#                                     protein_name, feature_type_name, desc, organism_name, 
#                                     ECO_str, source_str, ids_str, sequence])
            
#         tally += len(data["results"])
#         total = response.headers.get("x-total-results", "?")
#         logging.info("Retrieved %d of %s total results", tally, total)

#         url = response.links.get("next", {}).get("url")
#         if tally >= int(total):
#             break

#     # adds headers to dataframe
#     dataframe = pd.DataFrame(dataframe, columns = ["Unique_ID", "Accession", "Position", "Length", "Entry_type",
#                                                     "Protein", "Feature_type", "Description", "Organism",
#                                                     "ECO_codes", "Sources", "Source_ids", "Sequence"])

#     logging.info(f"Entry download is done. {len(dataframe)} records stored for further deduplication.")

#     ### for development purposes - can remove later
#     # stores all records
#     dataframe.to_csv(f"{data_dir}/all_records-{setup.version}.csv", index=False)

#     # deduplicates based on Accession + position + description and combines source info
#     collapsed = deduplicate_records(dataframe)

#     logging.info(f"Deduplication is done. {len(collapsed)} records stored for further processing.")

#     ### for development purposes - can remove later
#     # stores deduplicated records
#     collapsed.to_csv(f"{data_dir}/deduplicated_records-{setup.version}.csv", index=False)

#     return collapsed
