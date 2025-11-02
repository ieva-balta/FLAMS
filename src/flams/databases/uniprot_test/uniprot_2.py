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
MODIFICATIONS = {
    "acetylation": [r"acetyl\w*"],
    "adp-ribosylation": [r"ADP[-\s]?ribosyl[\w-]*"],
    "adp-riboxanation": [r"ADP[-\s]?ribox[\w-]*"], # new type
    "amidation": [r"amid[\w-]*"],
    "ampylation" : [r"AMP[\w-]*"],
    "benzoylation": [r"benzoyl[\w-]*"],
    "beta-hydroxybutyrylation": [r"(?:β|beta)[-\s]?hydroxybutyryl[\w-]*"],
    "biotinylation": [r"biotinyl[\w-]*"],
    "blocked_amino_end": [r"blocked amino end"],
    "bromination": [r"bromo[\w-]*"], # new type
    # butyrylation moved to the ends bc its more general
    # "carbamidation": "",
    "carboxyethylation": [r"carboxyethyl[\w-]*"],
    # carboxylation moved to broad
    # "carboxymethylation": "",
    "cholesterol_ester": [r"cholesterol(?: glycine)? ester"],
    "citrullination": [r"citrulline"],
    "crotonylation": [r"crotonyl[\w-]*"],
    "cyclopeptide": [r"cyclopeptide"],
    # cysteinylation moved to broad terms
    "c-linked_glycosylation": [r"\bc[-\s]?linked"],
    "deamidation": [r"deamidat\w*"],
    "deamination": [r"allysine"],
    "decanoylation": [r"decanoyl[\w-]*"],
    "decarboxylation": [r"decarboxylat\w*"],
    "dehydration": [r"dehydro[\w-]*"], # new type
    "dephosphorylation": [r"dephospho[\w-]*"],
    # "dietylphosphorylation": "",
    # leaving this for later:
    # "disulfide_bond": "ft_disulfid:*", # maybe ft_chain: disulfide bond as well? or smotheing liek it
    # d-enantiomers moved
    "d-glucuronoylation": [r"d[-\s]?glucuronoyl[\w-]*"],
    "farnesylation": [r"farnesyl[\w-]*"],
    # formation of an isopeptide moved to broad
    "formylation": [r"formyl[\w-]*"],
    "gamma-carboxyglutamic_acid": [r"(?:gamma|γ)-?carboxyglutamic acid"],
    #genarylation moved to broad
    "geranylgeranylation": [r"geranylgeranyl[\w-]*"], 
    "glutamylation": [r"glutamyl[\w-]*"], # new type
    "glutarylation": [r"glutaryl[\w-]*"],
    "glutathionylation": [r"glutathionyl(?: cysteine)?"],
    "glycation": [r"glycation"],
    "gmpylation": [r"GMP[\w-]*"], # new type
    "gpi-anchor": [r"gpi[-\s]?anchor"],
    "histidylation": [r"histidyl[\w-]*"], # new type
    # "hmgylation": "",
    "hydroxyceramide_ester": [r"hydroxyceramide(?: glutamate)? ester"],
    # hydroxylation moved to broad
    "hypusine": [r"hypusine"], # new type
    "imidazolation": [r"imidazol[\w-]*"], # new type
    "iodination": [r"iodo[\w-]*"],
    "lactoylation": [r"lactoyl[\w-]*"],
    # "lactylation": "",
    "lipoylation": [r"lipoyl[\w-]*"],
    "malonylation": [r"malonyl[\w-]*"],
    #methylation moved to broad
    # "mgcylation": "",
    # "mgylation": "",
    "myristoylation": [r"myristoyl[\w-]*"],
    "neddylation": [r"NEDD8"],
    # nitration moved to broad 
    "n-carbamoylation": [r"n[-\s]?carbamoyl[\w-]*"],
    "n-linked_glycosylation": [r"n[\s-]?(?:alpha|beta|α|β)?[\s-]?linked"],
    "n-palmitoylation": [r"n[\s-]?palmitoyl[\w-]*"],
    "octanoylation": [r"octanoyl[\w-]*"],
    # oxidation moved to broad terms
    "o-linked_glycosylation": [r"o[\s-]?(?:alpha|beta|α|β)?[\s-]?linked"],
    "o-palmitoleoylation": [r"o[\s-]?palmitoleoyl[\w-]*"],
    "o-palmitoylation": [r"o[\s-]?palmitoyl[\w-]*"], # overlap with palmitoylation
    # palmitoylation moved to broad
    "phosphatidylethanolamine_amidation": [r"phosphatidylethanolamine amidated"],
    # "phosphoglycerylation": "",
    #phosphorylation moved to broad
    "propionylation": [r"propionyl[\w-]*"],
    "pupylation": [r"Pup[\w-]*"],
    "pyridoxal_phosphate_addition": [r"pyridoxal phosphate"],
    "pyrrolidone_carboxylic_acid": [r"pyrrolidone carboxylic acid"],
    "pyrrolylation": [r"pyrrolyl[\w-]*"],
    "pyruvate": [r"pyruvic acid"],
    "serotonylation" : [r"seroton[\w-]*"],
    "stearoylation" : [r"stearoyl[\w-]*"],
    "succinylation": [r"succiny[\w-]*"],
    # sulfation moved to broad
    "sulfhydration" : [r"cysteine persulfide"],
    "sulfoxidation" : [r"sulfoxide"],
    "sumoylation": [r"SUMO"], # doesn't query "Modified residue (large scale data)" : "Sumoylated lysine"
    "s-archaeol" : [r"s[\s-]?archaeol[\w-]*"],
    "s-carbamoylation" : [r"s[\s-]?carbamoyl[\w-]*"],
    "s-cyanation" : [r"s[\s-]?cyano[\w-]*"],
    "s-cysteinylation" : [r"s[\s-]?cysteinyl[\w-]*"],
    "s-diacylglycerol" : [r"s[\s-]?diacylglycerol[\w-]*"],
    "s-linked_glycosylation": [r"s[\s-]?(?:alpha|beta|α|β)?[\s-]?linked"],
    "s-nitrosylation" : [r"s[\s-]?nitro[\w-]*"], # overlaps with nitro*
    "s-palmitoylation" : [r"s[\s-]?palmitoyl[\w-]*"],
    "thiocarboxylation" : [r"thioglycine"],
    "ubiquitination": [r"ubiquitin"], # same problem as sumoylation
    "umpylation" : [r"UMP[\w-]*"],
    "2-hydroxyisobutyrylation": [r"2[-\s]?hydroxyisobutyryl[\w-]*"], # overlap with butyryl*

    # broader terms at the end:
    "butyrylation": [r"butyryl[\w-]*"], # overlaps with beta-hydroxybutyrylation and 2-hydroxyisobutyrylation
    "cysteinylation": [r"cysteinyl[\w-]*"], # new type, overlap with s-cysteinylation etc
    "carboxylation": [r"carboxy[\w-]*"], #this might be a broad one, carboxylic acid and carboxylysine and carboxyglutamate
    # "d-enantiomers": "ft_mod_res:d-", # new type, is applicable? (reported as modified residues on uniprot)
    "formation_of_an_isopeptide_bond": [r"isopeptide"], # lots of overlap with the sumo, ubiq, nedd etc.
    "geranylation": [r"geranyl[\w-]*"], # new type
    "hydroxylation": [r"hydroxy[\w-]*"], # doesn't catch all dihydroxy
    "methylation": [r"methyl[\w-]*"],
    "oxidation": [r"oxo[\w-]*", r"sulfenic acid"], #Tryptophylquinone? Methionine sulfoxide? Cysteine sulfenic acid?
    "palmitoylation": [r"palmitoyl[\w-]*"], # new type
    "phosphorylation": [r"phospho[\w-]*"],
    "nitration": [r"nitro[\w-]*"],
    "sulfation" : [r"sulfo[\w-]*"],
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


# to not run the download all the time
df = pd.read_csv("/data/leuven/368/vsc36826/flams/flams2/FLAMS/src/flams/databases/uniprot_test/records.csv", header=0)

def get_fasta(PTM_modification_dict, data_dir):
    """
    make fasta files
    log unclassified records
    """

    # switch when integrated
#    unclassified = sort_uniprot_records(get_uniprot_records(), PTM_modification_dict, data_dir)
    unclassified = sort_uniprot_records(df, PTM_modification_dict, data_dir)

    #special case disulfide bonds: could move to the download bit before the dataframe is returned
    disulfide_bonds = unclassified[
        (unclassified["Feature_type"] == "Chain") &
        (unclassified["Description"] == "nan")
    ]
    disulfide_bonds["Description"] = "Disulfide__bond"
    fasta_records_disulfide_bonds = df_to_fasta(disulfide_bonds)
    with open(f"{data_dir}/disulfide_bonds-{str(version)}.fasta", "w", encoding="UTF-8") as out:
        SeqIO.write(fasta_records_disulfide_bonds, out, "fasta")
    logging.info(f"Fasta file with {len(disulfide_bonds)} records for modifictation disulfide bonds created and stored at {data_dir}.")

    # remove them from the original
    unclassified = unclassified[~(
        (unclassified["Feature_type"] == "Chain") &
        (unclassified["Description"] == "nan")
    )]

    # remove the ones with removed and interchain description:
    unclassified = unclassified[~(
        (unclassified["Description"] == "Interchain") &
        (unclassified["Description"] == "Removed")
    )]

    fasta_records_unclassified = df_to_fasta(unclassified)

    with open(f"{data_dir}/unclassified-{str(version)}.fasta", "w", encoding="UTF-8") as out:
        SeqIO.write(fasta_records_unclassified, out, "fasta")
    logging.info(f"Fasta file for unclassified entries was created and stored at {data_dir}/unclassified-{str(version)}.fasta.")
    logging.info(f"There are {len(unclassified)} unclassified entries.")


def sort_uniprot_records(uniprot_records, PTM_modification_dict, data_dir):
    """
    sort through uniprot records to group them by PTM type
    outputs fasta per modification type
    return unclassified records
    """
    # copy input dataframe
    remaining_records = uniprot_records.copy()

    for modification, regex in MODIFICATIONS.items():
        df_mod_type = find_modification_type(remaining_records, regex)

        #skip empty mathces
        if df_mod_type.empty:
            logging.info(f"No entries for modification type {modification} were found.")
            continue

        #removes classified records
        remaining_records = remaining_records.drop(df_mod_type.index)

        fasta_records = df_to_fasta(df_mod_type)
        
        with open(f"{data_dir}/{modification}-{str(version)}.fasta", "w", encoding="UTF-8") as out:
            SeqIO.write(fasta_records, out, "fasta")
        
        logging.info(f"Fasta file with {len(df_mod_type)} records for modifictation {modification} created and stored at {data_dir}.")
    
    return remaining_records
    

def find_modification_type(uniprot_records, PTM_regex_list):
    """
    returns dataframe with only records matching that PTM type
    """
    # Combine all regex patterns into one (using OR)
    combined_pattern = "|".join(f"(?:{p})" for p in PTM_regex_list)
    regex = re.compile(combined_pattern, flags=re.IGNORECASE)

    # Use vectorized str.contains() to filter rows
    mask = uniprot_records["Description"].str.contains(regex, na=False)
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
    
def get_uniprot_records():
    """
    fetches all uniprot ptm records
    outputs a dataframe of all records + info
    """

    logging.info(f"Fetching entries from UniProt, please wait.")

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url = f"{base_url}?query=existence:1&format=json"

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
            name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")
            entry_type = entry.get("entryType", "")

            # to check for duplicated entries down the line
            ptm_sites = []

            # find PTM sites
            for feature in entry.get("features", []):
                feature_type = feature.get("type")
                if feature.get("type") not in ["Modified residue",
                            "Chain",
                            "Modified residue (large scale data)", 
                            "Lipidation", 
                            "Glycosylation",
                            "Disulfide bond",
                            "Cross-link",
                            "Initiator methionine"]: #prob remove initiator methionine
                    continue

                desc = feature.get("description", "N/A") 
                
                # remove the notes from the description
                if ";" in desc:
                    desc = desc.split(";")[0].strip()

                #get psoition
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                
                #check duplicate ptm site entries
                ptm_desc = f"{pos}|{desc}"
                if ptm_desc in ptm_sites:
                    continue

                ptm_sites.append(ptm_desc)

                # get evidence ECO|source|ids
                ECOs = []
                sources = []
                ids = []
                for ev in feature.get("evidences", []):
                    eco = ev.get("evidenceCode", "")
                    # skips if the evidnece is from similar protein
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

    logging.info(f"Entry download is done. {len(dataframe)} records stored for further processing.")

    #remove later
    dataframe.to_csv("records.csv", index=False)

    return dataframe

# run the thing
get_fasta(MODIFICATIONS, "/data/leuven/368/vsc36826/flams/flams2/dbs")