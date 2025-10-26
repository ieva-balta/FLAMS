#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" uniprot
functions to fetch and parse uniprot entries
based on scripts for cplm4 and dbptm
"""
# to run example on bash:
# python -c "import uniprot; uniprot.get_fasta_rest('acetyl', '/data/leuven/368/vsc36826/flams/flams2/dbs/acetyl_rest.fasta')"

import sys
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#fuzzywuzzy to add to conda download?
import requests
from fuzzywuzzy import process
import logging

#might remove
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def get_fasta_rest(descriptor, location):
    """
    will download uniprot data for each descriptor (modification) seperately
    and save it in a fasta format
    return the multi fasta file for the specified decriptor
    """
    
    # split if multiple keywords used
    # remove the feature type
    if "OR" in descriptor:
        descrip = descriptor.replace("(", "").replace(")", "")
        descrip = descrip.split(" OR ")
        descriptors = [keyword.split(":")[1].replace("*", "").replace("-", " ") for keyword in descrip]
    else:
        descriptors = [descriptor.split(":")[1].replace("*", "").replace("-", " ")]

    descriptor_text = ", ".join(descriptors)
    logging.info(f"Downloading UniProt {descriptor_text} Database, please wait.")

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url = f"{base_url}?query={descriptor} AND existence:1&format=json"

    tally = 0
    fasta_records = []

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

            # to check for duplication down the line
            ptm_sites = []

            # find PTM sites
            for feature in entry.get("features", []):
                if feature.get("type") not in ["Modified residue", 
                            "Modified residue (large scale data)", 
                            "Lipidation", 
                            "Glycosylation",
                            "Disulfide bond",
                            "Cross-link",
                            "Initiator methionine"]:
                    continue

                desc = feature.get("description", "N/A") 
                
                # remove the notes from the description
                if ";" in desc:
                    desc = desc.split(";")[0].strip()

                # special case - disulfide bonds
                if descriptor == "ft_disulfid:*" and (feature.get("type") != "Disulfide bond" or "interchain" in desc.lower()):
                    continue

                #try to match them ignoring the extra characthers like * and -
                # extra step to continue the for feature loop if there ar eno matches
                skip_outer_for_loop = True # default True
                if descriptor != "ft_disulfid:*":
                    for d in descriptors:
                        if d.lower() in desc.lower().replace("-"," "):
                            skip_outer_for_loop = False # flip to False if there is a match
                            break
                if skip_outer_for_loop and descriptor != "ft_disulfid:*": # skip if no match
                    continue

                # special case - disulfide bonds
                if descriptor == "ft_disulfid:*":
                    desc = "Disulfide bond"

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
                    if eco not in ["ECO:0000269", "ECO:0000314", "ECO:0007744", "ECO:0007829",
                                   "ECO:0000312", "ECO:0000313"]:
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

                #write the fasta record for each PTM site
                entry_type = entry_type.replace(" ", "__")
                protein_name = f"{name}".replace(" ","__")
                organism_name = f"{organism}".replace(" ","__")
                desc = f"{desc}".replace(" ","__")
                seq = Seq(sequence)
                length = len(seq)
                id = f"{accession}|{pos}|{length}|{entry_type}"
                rec = SeqRecord(
                    seq,
                    id=id,
                    description=f"{protein_name}|{desc}|{organism_name} [{ECO_str}|{source_str}|{ids_str}]",
                )

                #append
                fasta_records.append(rec)
            
        tally += len(data["results"])
        total = response.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results for descriptor %s", tally, total, descriptor)

        url = response.links.get("next", {}).get("url")
        if tally >= int(total):
            break

    with open(location, "a", encoding="UTF-8") as out:
        SeqIO.write(fasta_records, out, "fasta")

    logging.info(f"Converted and stored UniProt {descriptor} Database entries as FASTA entries for the local {descriptor} BLAST database format.")
    logging.info(f"Total of {len(fasta_records)} fasta records stored for descriptor {descriptor}.")