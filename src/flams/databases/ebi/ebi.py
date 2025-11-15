import sys
import json
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import requests
import logging

import time

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


#valid eco codes
valid_ECO_codes = [
    "ECO:0000269", #(experimental evidence used in manual assertion)​
    "ECO:0000314", #(direct assay evidence used in manual assertion)​
    "ECO:0007744", #(combinatorial computational and experimental evidence used in manual assertion)​
    "ECO:0007829", #(combinatorial computational and experimental evidence used in automatic assertion)​
]

def get_uniprot_records(data_dir):
    """
    fetches all uniprot ptm records
    outputs a dataframe of all records + info

    This function fetches all relevant uniprot entries (existence:1 and valid_eco_codes) and stores them in a pandas dataframe.

    Parameters
    ----------
    data_dir: str
        Location for storing data

    """
    logging.info(f"Fetching entries from UniProt, please wait.")


    base_url_rest = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url_rest = f"{base_url_rest}?query=existence:1&format=json&size=500"

    # for pagination purposes
    tally = 0

    # to store data
    dataframe = []

    # to keep track of modified residues all together
    tally_per_feature = 0

    while url_rest:
        response = requests.get(url_rest, headers=headers)
        response.raise_for_status()
        data = response.json()

        # parse through results
        for entry in data.get("results", []):

            start_time = time.time()

            accession = entry.get("primaryAccession")
            sequence = entry.get("sequence", {}).get("value", "")
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

                # gets evidence ECO|source|ids
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

                # removes white space
                feature_type_name = f"{feature_type}".replace(" ", "__")
                desc = f"{desc}".replace(" ","__")
                
                # adds ID to differentiate
                unique_id = f"{accession}_{tally_per_feature}"
                tally_per_feature += 1

                # appends to dataframe
                # add with start position
                dataframe.append([unique_id, accession, pos_start, length, entry_type, 
                                protein_name, feature_type_name, desc, organism_name, 
                                ECO_str, source_str, ids_str, 
                                "", "", "", "", "",
                                sequence])

                # adds ID to differentiate
                unique_id = f"{accession}_{tally_per_feature}"
                tally_per_feature += 1
                
                # adds the end position
                dataframe.append([unique_id, accession, pos_end, length, entry_type, 
                                protein_name, feature_type_name, desc, organism_name, 
                                ECO_str, source_str, ids_str, 
                                "", "", "", "", "",
                                sequence])
            # find PTMS in EBI API            
            url_ebi = f"https://www.ebi.ac.uk/proteins/api/proteomics/ptm/{accession}"

        
            try:
                response_ebi = requests.get(url_ebi, headers=headers)
                response_ebi.raise_for_status()
            except requests.exceptions.HTTPError as e:
                # logging.info(f"Entry with accession {accession} not found on EBI API.")
                continue

            data_ebi = response_ebi.json()


            for feature in data_ebi.get("features", []):
                feature_type = feature.get("type", "")
                if feature_type not in ["PROTEOMICS_PTM"]:
                    continue

                # peptide start position
                pep_start = int(feature.get("begin", 0))

                #eco codes
                ecos = []
                for evidence in feature.get("evidences", []):
                    eco = evidence.get("code","")
                    if eco not in valid_ECO_codes:
                        continue
                    ecos.append(eco)
                ecos_str = ";".join(ecos)

                # PTM description
                for ptm in feature.get("ptms", []):
                
                    #description
                    desc = ptm.get("name", "").replace(" ", "__")

                    # get modified amino acid position 
                    pos_in_peptide = int(ptm.get("position", 0))
                    position = pep_start + pos_in_peptide - 1

                    sources = ";".join(ptm.get("sources", []))

                    #ptmxchange ids and sources and confidence scores
                    x_ids = []
                    pubmed_ids = []
                    conf_scores = []
                    for reference in ptm.get("dbReferences", []):
                        x_ids.append(reference.get("id", ""))
                        pubmed_ids.append(reference.get("properties", {}).get("Pubmed ID", ""))
                        conf_scores.append(reference.get("properties", {}).get("Confidence score", ""))
                    x_ids_str = ";".join(x_ids)
                    pubmed_ids_str = ";".join(pubmed_ids)
                    conf_scores_str = ";".join(conf_scores)

                    # adds ID to differentiate
                    unique_id = f"{accession}_{tally_per_feature}"
                    tally_per_feature += 1
                
                    # adds the end position
                    dataframe.append([unique_id, accession, position, length, entry_type, 
                            protein_name, feature_type, desc, organism_name,
                            "", "", "", 
                            ecos_str, sources, x_ids_str, pubmed_ids_str, conf_scores_str,
                            sequence])
            end_time = time.time()
            elapsed = end_time - start_time

            logging.info(f"Time elapsed for {accession}: {elapsed:.2f} seconds.")

        tally += len(data["results"])
        total = response.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)

        url_rest = response.links.get("next", {}).get("url")
        if tally >= int(total):
            break

    # adds headers to dataframe
    dataframe = pd.DataFrame(dataframe, columns = ["Unique_ID", "Accession", "Position", "Length", "Entry_type",
                                                    "Protein", "Feature_type", "Description", "Organism",
                                                    "ECO_codes", "Sources", "Source_ids", # Uniprot
                                                    "ECO_codes_LSD", "Sources_LSD", "LSD_IDs", "PubMed_LSD_IDs", "Confidence_scores_LSD", #PTMxChange or Pride
                                                    "Sequence"])

    logging.info(f"Entry download is done. {len(dataframe)} records stored for further deduplication.")

    ### for development purposes - can remove later
    # stores all records
    dataframe.to_csv(f"{data_dir}/O14929.csv", index=False)

    # # deduplicates based on Accession + position + description and combines source info
    # collapsed = deduplicate_records(dataframe)

    # logging.info(f"Deduplication is done. {len(collapsed)} records stored for further processing.")

    # ### for development purposes - can remove later
    # # stores deduplicated records
    # collapsed.to_csv(f"{data_dir}/deduplicated_records-{setup.version}.csv", index=False)

    # return collapsed

get_uniprot_records("/data/leuven/368/vsc36826/scratch/IBP/ebi")

