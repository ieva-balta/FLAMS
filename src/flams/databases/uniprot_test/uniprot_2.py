#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
will fetch all ptm records in one go and filter after
"""

import sys
import json
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#fuzzywuzzy to add to conda download?
import requests
from fuzzywuzzy import process
import logging

# regular expression re
import re

#might remove
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# dictionary of PTM types and a list of RegEx's to find them
# ALLOWS DUPLICATION
MODIFICATIONS = {
    "acetylation": [r"[-\w_]*acetyl[-\w_]*"],
    "adp-ribosylation": [r"ADP[-_\s]*ribosyl[-\w_]*"],
    "adp-riboxanation": [r"ADP[-_\s]*ribox[-\w_]*"], # new type
    "amidation": [r"[-\w_]*(?<!de)amid[-\w_]*"],
    "ampylation" : [r"AMP[-\w_]*", r"[-\w_]*adenylate[-\w_]*"],
    "benzoylation": [r"[-\w_]*benzoyl[-\w_]*"],
    "beta-hydroxybutyrylation": [r"(?:β|beta)[-_\s]*hydroxybutyryl[-\w_]*"], # overlaps with butyrylation and hydroxylation
    "biotinylation": [r"[-\w_]*biotinyl[-\w_]*"],
    "blocked_amino_end": [r"blocked__(?:amino|carboxyl)__end"],
    "bromination": [r"[-\w_]*bromo[-\w_]*"], # new type
    "butyrylation": [r"[_\w-]*butyryl[_\w-]*"], # some entries overlap with beta-hydroxybutyrylation and 2-hydroxyisobutyrylation
    # "carbamidation": "",
    "carboxyethylation": [r"[_\w-]*carboxyethyl[_\w-]*"],
    "carboxylation": [r"[_\w-]*(?<!de)carboxy[_\w-]*"], # overlaps
    # "carboxymethylation": "",
    "cholesterol_ester": [r"[-\w_]*cholesterol(?:__glycine)?__ester"],
    "citrullination": [r"[-\w_]*citrulline"],
    "crotonylation": [r"[-\w_]*crotonyl[-\w_]*"],
    "cyclopeptide": [r"[-\w_]*cyclo(?:peptide)?[-\w_]*"],
    "cysteinylation": [r"[_\w-]*cysteinyl[_\w-]*"], # overlap with s-cysteinylation etc
    "c-linked_glycosylation": [r"c[-_\s]*linked[-\w_]*"],
    "deamidation": [r"[-\w_]*deamidat[-\w_]*"],
    "deamination": [r"[-\w_]*allysine[-\w_]*"],
    "decanoylation": [r"[-\w_]*decanoyl[-\w_]*"],
    "decarboxylation": [r"[-\w_]*decarboxylat[-\w_]*"],
    "dehydration": [r"[-\w_]*dehydro[-\w_]*"], # new type
    "dephosphorylation": [r"[-\w_]*dephospho[-\w_]*"],
    # "dietylphosphorylation": "",
    "disulfide_bond": [r"[-\w_]*disulfide[-\w_]*"],
    "d-glucuronoylation": [r"d[-_\s]*glucuronoyl[-\w_]*"],
    "FADylation" : [r"FAD[-\w_]*"], #new type
    "farnesylation": [r"[-\w_]*farnesyl[-\w_]*"],
    "formation_of_an_isopeptide_bond": [r"[_\w-]*isopeptide[_\w-]*"], # lots of overlap with the sumo, ubiq, nedd etc.
    "formylation": [r"[-\w_]*formyl[-\w_]*"],
    "gamma-carboxyglutamic_acid": [r"(?:gamma|γ)-?carboxyglutamic__acid[-\w_]*"],
    #genarylation moved to broad
    "geranylation": [r"[_\w-]*geranyl[_\w-]*"], # new type, overlap with genarylgenarylation
    "geranylgeranylation": [r"[-\w_]*geranylgeranyl[-\w_]*"], 
    "glutamylation": [r"[-\w_]*glutamyl[-\w_]*"], # new type, some overlap with serotonin etc
    "glutarylation": [r"[-\w_]*glutaryl[-\w_]*"],
    "glutathionylation": [r"[-\w_]*glutathionyl(?:__cysteine)?[-\w_]*"],
    "glycation": [r"[-\w_]*glycation[-\w_]*"],
    "gmpylation": [r"[-\w_]*GMP[-\w_]*"], # new type
    "gpi-anchor": [r"[-\w_]*gpi[-\s]?(?:like__)anchor[-\w_]*"], # adds mimick "gpi-like anchor"
    "histidylation": [r"[-\w_]*histidyl[-\w_]*"], # new type
    # "hmgylation": "",
    "hydroxyceramide_ester": [r"[-\w_]*hydroxyceramide(?:__glutamate)?__ester[-\w_]*"],
    "hydroxylation": [r"[_\w-]*(?<!de)hydroxy[_\w-]*"],
    "hypusine": [r"[-\w_]*hypusine[-\w_]*"], # new type
    "imidazolation": [r"[-\w_]*imidazol[-\w_]*"], # new type
    "iodination": [r"[-\w_]*iodo[_\w-]*", r"[_\w-]*thyroxine[_\w-]*"],
    "isomerization": [r"[dl]-[_\w-]*"], # new type
    "lactoylation": [r"[_\w-]*lactoyl[_\w-]*"],
    # "lactylation": "",
    "lipoylation": [r"[_\w-]*lipoyl[_\w-]*"],
    "malonylation": [r"[_\w-]*malonyl[_\w-]*"],
    "methylation": [r"[_\w-]*methyl[_\w-]*"],
    # "mgcylation": "",
    # "mgylation": "",
    "myristoylation": [r"[_\w-]*myristoyl[_\w-]*"],
    "neddylation": [r"[_\w-]*NEDD8[_\w-]*"],
    "nitration": [r"[_\w-]*nitro[_\w-]*", r"[_\w-]*nitrated[_\w-]*"],
    "n-carbamoylation": [r"n[-_\s]*carbamoyl[_\w-]*"],
    "n-linked_glycosylation": [r"n[-_\s]*(?:alpha|beta|α|β)?[\s-]?linked[_\w-]*"],
    "n-palmitoylation": [r"n[-_\s]*palmitoyl[_\w-]*"],
    "octanoylation": [r"[_\w-]*octanoyl[_\w-]*"],
    "oxidation": [r"[_\w-]*oxo[_\w-]*", r"[_\w-]*sulf[ei]nic__acid[_\w-]*"], #Tryptophylquinone? Methionine sulfoxide? Cysteine sulfenic acid?
    "o-linked_glycosylation": [r"o[-_\s]*(?:alpha|beta|α|β)?[\s-]?linked[_\w-]*"],
    "o-palmitoleoylation": [r"o[-_\s]*palmitoleoyl[_\w-]*"],
    "o-palmitoylation": [r"o[-_\s]*palmitoyl[_\w-]*"], # overlap with palmitoylation
    "palmitoylation": [r"[_\w-]*palmitoyl[_\w-]*"], # new type
    "phosphatidylethanolamine_amidation": [r"[_\w-]*phosphatidylethanolamine__amidated[_\w-]*"],
    "phosphoglycerylation": [r"[_\w-]*glycerophospho[_\w-]*"], # overlaps
    "phosphorylation": [r"[_\w-]*(?<!de)phospho[_\w-]*", r"[_\w-]*aspartylphosphate[_\w-]*"],
    "prenylation" : [r"[_\w-]*prenyl[_\w-]*"],
    "propionylation": [r"[_\w-]*propionyl[_\w-]*"],
    "pupylation": [r"[_\w-]*pup[_\w-]*"],
    "pyridoxal_phosphate_addition": [r"[_\w-]*pyridoxal__phosphate[_\w-]*"],
    "pyrrolidone_carboxylic_acid": [r"[_\w-]*pyrrolidone__carboxylic__acid[_\w-]*"], # overlap
    "pyrrolylation": [r"[_\w-]*pyrrolyl[_\w-]*"],
    "pyruvate": [r"[_\w-]*pyruvic__acid[_\w-]*"], # overlap
    "quinones" : [r"[_\w-]*quinone[_\w-]*"],
    "serotonylation" : [r"[_\w-]*seroton[_\w-]*"], # overlap
    "stearoylation" : [r"[_\w-]*stearoyl[_\w-]*"],
    "succinylation": [r"[_\w-]*succin[iy][_\w-]*"],
    "sulfation" : [r"[_\w-]*sulfo[_\w-]*"],
    "sulfhydration" : [r"[_\w-]*cysteine__persulfide[_\w-]*"],
    "sulfilimine_crosslink" : [r"[_\w-]*sulfilimine[_\w-]*"], # new type
    "sulfoxidation" : [r"[_\w-]*sulfoxide[_\w-]*"],
    "sumoylation": [r"[_\w-]*SUMO[_\w-]*"], 
    "s-archaeol" : [r"s[_\s-]?archaeol[_\w-]*"],
    "s-carbamoylation" : [r"s[-_\s]*carbamoyl[_\w-]*"],
    "s-cyanation" : [r"s[-_\s]*cyano[_\w-]*"],
    "s-cysteinylation" : [r"s[-_\s]*cysteinyl[_\w-]*"],
    "s-diacylglycerol" : [r"s[-_\s]*diacylglycerol[_\w-]*"],
    "s-linked_glycosylation": [r"s[-_\s]*(?:alpha|beta|α|β)?[\s-]?linked[_\w-]*"],
    "s-nitrosylation" : [r"s[-_\s]*nitro[_\w-]*"], # overlaps with nitro*
    "s-palmitoylation" : [r"s[-_\s]*palmitoyl[_\w-]*"],
    "thiocarboxylation" : [r"[_\w-]*thioglycine[_\w-]*"],
    "thioester_crosslink" : [r"[_\w-]*thioester[_\w-]*"],
    "ubiquitination": [r"[_\w-]*ubiquitin[_\w-]*"],
    "umpylation" : [r"UMP[_\w-]*"],
    "2-hydroxyisobutyrylation": [r"2[-_\s]*hydroxyisobutyryl[_\w-]*"], # overlap with butyrylation and hydroxylation
    }

#valid eco codes
valid_ECO_codes = [
    "ECO:0000269", #(experimental evidence used in manual assertion)​
    "ECO:0000314", #(direct assay evidence used in manual assertion)​
    "ECO:0007744", #(combinatorial computational and experimental evidence used in manual assertion)​
    "ECO:0007829", #(combinatorial computational and experimental evidence used in automatic assertion)​
]

#database version
version = 2.0


# # to not run the download all the time
# df = pd.read_csv("/data/leuven/368/vsc36826/flams/flams2/FLAMS/src/flams/databases/uniprot_test/records.csv", header=0)

def get_fasta(PTM_modification_dict, data_dir):
    """
    make fasta files
    log unclassified records
    """

    # switch when integrated
    # all_records = df
    all_records = get_uniprot_records()

    # classifies records into PTM types and stores fasta files per modification type
    classified = sort_uniprot_records(all_records, PTM_modification_dict, data_dir)

    # find unclassified
    unclassified = pd.concat([all_records, classified]).drop_duplicates(keep=False)

    #remove later
    classified.to_csv("classified.csv", index=False)
    unclassified.to_csv("unclassified.csv", index=False)

    #output unclassified in fasta format
    fasta_records_unclassified = df_to_fasta(unclassified)

    with open(f"{data_dir}/unclassified-{str(version)}.fasta", "w", encoding="UTF-8") as out:
        SeqIO.write(fasta_records_unclassified, out, "fasta")

    logging.info(f"Fasta file for {len(unclassified)} unclassified entries was created and stored at {data_dir}/unclassified-{str(version)}.fasta.")


def sort_uniprot_records(uniprot_records, PTM_modification_dict, data_dir):
    """
    sort through uniprot records to group them by PTM type
    outputs fasta per modification type
    return unclassified records
    """
    # copy input dataframe
    classified_records = pd.DataFrame(columns = ["Accession", "Position", "Length", "Entry_type",
                                                    "Protein", "Feature_type", "Description", "Organism",
                                                    "ECO_codes", "Sources", "Source_ids", "Sequence"])

    for modification, regex in MODIFICATIONS.items():
        df_mod_type = find_modification_type(uniprot_records, regex)

        #skip empty mathces
        if df_mod_type.empty:
            logging.info(f"No entries for modification type {modification} were found.")
            continue

        #records classified records
        classified_records = pd.concat([classified_records, df_mod_type], ignore_index=True)

        fasta_records = df_to_fasta(df_mod_type)

        #remove later
        df_mod_type.to_csv(f"modification_csvs/{modification}.csv", index=False)
        
        with open(f"{data_dir}/{modification}-{str(version)}.fasta", "w", encoding="UTF-8") as out:
            SeqIO.write(fasta_records, out, "fasta")
        
        logging.info(f"Fasta file with {len(df_mod_type)} records for modifictation {modification} created and stored at {data_dir}.")
    
    return classified_records
    

def find_modification_type(uniprot_records, PTM_regex_list):
    """
    returns dataframe with only records matching that PTM type
    """
    # Combine all regex patterns into one (using OR)
    combined_pattern = "|".join(f"(?:{p.lower()})" for p in PTM_regex_list)

    # Use vectorized str.contains() to filter rows
    description = uniprot_records["Description"].astype(str).str.lower()
    mask = description.str.contains(combined_pattern, na=False, regex=True, flags=re.IGNORECASE)
    df_filtered = uniprot_records.loc[mask].copy()

    return df_filtered

def df_to_fasta(PTM_type_df):
    """
    formats dataframe into a fasta records dataframe
    """
    fasta_records = []
    for idx, row in PTM_type_df.iterrows():
        seq = Seq(row["Sequence"])
        id = f'{row["Accession"]}|{str(row["Position"])}|{row["Length"]}|{row["Entry_type"]}'
        rec = SeqRecord(
                    seq,
                    id=id,
                    description=f'{row["Protein"]}|{row["Description"]}|{row["Organism"]} [{row["ECO_codes"]}|{row["Sources"]}|{row["Source_ids"]}]',
                )

        #append
        fasta_records.append(rec)

    return fasta_records

def deduplicate_records(uniprot_records):
    """
    takes dataframe of all uniprot PTM records
    deduplicates based on uniprot accesion + PTM description + position
    concatinates ECO codes and sources
    """

    # columns to group by
    group_cols = ["Accession", "Position", "Description"]

    collapsed = (
        uniprot_records.groupby(group_cols, as_index=False)
        .agg({
          "Length": "first",
          "Entry_type": "first",
          "Protein": "first",
          "Feature_type": "first",
          "Organism": "first",
          "Sequence": "first",
          # combine these fields by joining unique non-null values
          "ECO_codes": lambda x: ";".join(sorted(set(filter(None, x)))),
          "Sources": lambda x: ";".join(sorted(set(filter(None, x)))),
          "Source_ids": lambda x: ";".join(sorted(set(filter(None, x))))
        })
    )

    return collapsed
    
def get_uniprot_records():
    """
    fetches all uniprot ptm records
    outputs a dataframe of all records + info
    """

    logging.info(f"Fetching entries from UniProt, please wait.")

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url = f"{base_url}?query=existence:1&format=json&size=500"

    tally = 0
    dataframe = []

    while url:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        data = response.json()

       # parse through results
        for entry in data.get("results", []):
            accession = entry.get("primaryAccession")
            sequence = entry.get("sequence", {}).get("value", "")
            #switch to gene name instead of protein name
            # name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            name = entry.get("genes", {}).get("geneName", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")

            entry_type = entry.get("entryType", "")
            if "Swiss-Prot" in entry_type:
                entry_type == "Swiss-Prot"
            if "TrEMBL" in entry_type:
                entry_type == "TrEMBL"

            # # to check for duplicated entries at the same position down the line
            # ptm_sites = []

            # find PTM sites
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

                desc = feature.get("description", "")

                #special case disulfide bonds
                if feature_type == "Disulfide bond":
                    desc = "Disulfide bond"

                # remove the notes from the description
                if ";" in desc:
                    desc = desc.split(";")[0].strip()

                #special case "removed"
                if desc.lower() == "removed":
                    continue

                #special case - microbial infection
                if "microbial infection" in desc.lower():
                    desc = desc.replace("(Microbial infection) ", "")

                #get psoition
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                
                # #check duplicate ptm site entries
                # ptm_desc = f"{pos}|{desc}"
                # if ptm_desc in ptm_sites:
                #     continue
                # ptm_sites.append(ptm_desc)

                # get evidence ECO|source|ids
                ECOs = []
                sources = []
                ids = []
                for ev in feature.get("evidences", []):
                    eco = ev.get("evidenceCode", "")
                    # skips if the evidnece is from the valid eco codes
                    if eco not in valid_ECO_codes:
                        continue
                    source = ev.get("source", "")
                    source_id = str(ev.get("id", ""))
                    ECOs.append(eco)
                    sources.append(source)
                    ids.append(source_id)
                ECO_str = ";".join(ECOs) 
                # skips if there are no ECOs
                if not ECOs:
                    continue 
                source_str = ";".join(sources) if sources else ""
                ids_str = ";".join(ids) if ids else ""

                #remove white space
                entry_type = entry_type.replace(" ", "__")
                protein_name = f"{name}".replace(" ","__")
                organism_name = f"{organism}".replace(" ","__")
                feature_type_name = f"{feature_type}".replace(" ", "__")
                desc = f"{desc}".replace(" ","__")
                length = len(sequence)

                #append to dataframe
                dataframe.append([accession, pos, length, entry_type, 
                                    protein_name, feature_type_name, desc, organism_name, 
                                    ECO_str, source_str, ids_str, sequence])
            
        tally += len(data["results"])
        total = response.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)

        url = response.links.get("next", {}).get("url")
        if tally >= int(total):
            break

    dataframe = pd.DataFrame(dataframe, columns = ["Accession", "Position", "Length", "Entry_type",
                                                    "Protein", "Feature_type", "Description", "Organism",
                                                    "ECO_codes", "Sources", "Source_ids", "Sequence"])

    logging.info(f"Entry download is done. {len(dataframe)} records stored for further deduplication.")

    # deduplicates based on Accession + position + description and combines source info
    collapsed = deduplicate_records(dataframe)

    logging.info(f"Deduplication is done. {len(collapsed)} records stored for further processing.")

    #remove later
    collapsed.to_csv("records.csv", index=False)

    return collapsed

# run the thing
get_fasta(MODIFICATIONS, "/data/leuven/368/vsc36826/flams/flams2/dbs")