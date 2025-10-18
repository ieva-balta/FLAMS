#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" uniprot
functions to fetch and parse uniprot entries
based on scripts for cplm4 and dbptm
"""
# to run example on bash:
# python -c "import uniprot; uniprot.get_fasta_rest('acetyl', '/data/leuven/368/vsc36826/flams/flams2/dbs/acetyl_rest.fasta')"

from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##################
# SPARQL version
##################


def get_fasta_sparql(descriptor, location):
    """
    will download uniprot data for each descriptor (modification) seperately
    and save it in a fasta format
    return the multi fasta file for the specified decriptor
    """
    
    with open(location, "a", encoding="UTF-8") as out:
        SeqIO.write(_convert_uniprot_to_fasta_sparql(_fetch_uniprot_sparql(descriptor)), out, "fasta")


### need to figure out what to do with different PTM names
def _fetch_uniprot_sparql(descriptor):
    """
    will run the sparql query
    putput is in json format

    Parameters
    ----------
    descriptor: str
        Description of a specific modification
    """

    endpoint = "https://sparql.uniprot.org/sparql"
    sparql = SPARQLWrapper(endpoint) # UniProt SPARQL endpoint
    sparql.setReturnFormat(JSON) # return in JSON format because it's easier to parse
    sparql.setTimeout(300)


    query = f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?accession ?position ?ptm_note ?aaSequence
    WHERE {{
        ?protein a up:Protein . # select proteins
        ?protein up:reviewed true . # only reviewed swiss-Prot entries
        ?protein up:sequence ?sequence . # get the sequence
        ?protein up:annotation ?ann . # get annotations
        ?ann a up:Modified_Residue_Annotation . # filter for modified residues
        ?ann rdfs:comment ?ptm_note . # get the PTM description
        ?ann up:range ?range . # get the range of the modification
        ?range faldo:begin ?loc . # get the start location
        ?loc faldo:position ?position . # get the position
        ?sequence rdf:value ?aaSequence . # get the amino acid sequence
        BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession) # extract accession from URI
        FILTER(CONTAINS(LCASE(?ptm_note), "{descriptor.lower()}")) # filter by PTM type
  
     }}
    ORDER BY ?accession ?position
  """

    # Execute the query
    try:
        sparql.setQuery(query)
        results = sparql.query().convert()
        #logging.info(f"Downloading UniProt {descriptor} Database, please wait.") # not sure how to add the size
    except Exception as e:
        print(f"Query failed with exception: {e}")
        return []

 
    # Extract and return results
    ptm_sites = []
    for row in results["results"]["bindings"]:
        sequence = row["aaSequence"]["value"]
        ptm_sites.append({
        "accession": row["accession"]["value"],
        "position": row["position"]["value"],
        "ptm_note": row["ptm_note"]["value"],
        "sequence": sequence
       })

    return ptm_sites

def _convert_uniprot_to_fasta_sparql(uniprot):
    """
    will help convert the query output to a fasta format
    returns sequence and records
    """
    fasta_records = []
    for site in uniprot:
        seq = Seq(site['sequence'])
        length = len(seq)
        id = f"{site['accession']}|{site['position']}|{length}|UniProt"
        rec = SeqRecord(
            seq, 
            id=id,
            description=f"{site['ptm_note']}")
        fasta_records.append(rec)
    return fasta_records


##################
# REST version
##################

#fuzzywuzzy to add to conda download?
import requests
from fuzzywuzzy import process
import logging

def get_fasta_rest(descriptor, location):
    """
    will download uniprot data for each descriptor (modification) seperately
    and save it in a fasta format
    return the multi fasta file for the specified decriptor
    """

#     with open(location, "a", encoding="UTF-8") as out:
#         ptm_data, protein_data = _fetch_uniprot_rest(descriptor)
#         SeqIO.write(_convert_uniprot_to_fasta_rest(ptm_data, protein_data), out, "fasta")


# def _fetch_uniprot_rest(descriptor):
#     # url
#     base_url = "https://rest.uniprot.org/uniprotkb/search"
#     headers = {"Accept": "application/json"}
#     url = f"{base_url}?query=reviewed:true AND ft_mod_res:{descriptor} AND existence:1&format=json"

#     tally = 0
#     protein_data = {}
#     ptm_data = []
#     while url:
#         response = requests.get(url, headers=headers)
#         response.raise_for_status()
#         data = response.json()

#        # parse through results
#         for entry in data.get("results", []):
#             accession = entry.get("primaryAccession")

#             # find PTM site
#             for feature in entry.get("features", []):
#                 if feature.get("type") != "Modified residue":
#                     continue
#                 desc = feature.get("description", "N/A")
#                 if descriptor.lower() not in desc.lower():
#                     continue
#                 pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
#                 evidences = [ev.get("evidenceCode") for ev in feature.get("evidences", [])]
#                 eco_code = evidences[0] if evidences else "N/A"

#                 ptm_data.append({
#                     "accession": accession,
#                     "position": pos,
#                     "ptm_note": desc,
#                     "evidence": eco_code
#                 })

#             if accession not in protein_data.keys():
#                 sequence = entry.get("sequence", {}).get("value", "")
#                 name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
#                 organism = entry.get("organism", {}).get("scientificName", "")
#                 protein_data[accession] = {
#                     "name": name,
#                     "organism": organism,
#                     "sequence": sequence,
#                 }

            
#         tally += len(data["results"])
#         total = response.headers.get("x-total-results", "?")

#         url = response.links.get("next", {}).get("url")
#         if tally >= int(total):
#             break

#     return ptm_data, protein_data

# def _convert_uniprot_to_fasta_rest(ptm_data, protein_data):
#     """
#     takes the ptm list and protein data and converts to fasta records
#     """
#     fasta_records = []
#     for ptm in ptm_data:
#         accession = ptm["accession"]
#         protein = protein_data.get(accession, {})
#         protein_name = f"{protein.get('name', '')}".replace(" ","__")
#         organism_name = f"{protein.get('organism', '')}".replace(" ","__")
#         seq = Seq(protein.get("sequence", ""))
#         length = len(seq)
#         id = f"{accession}|{ptm['position']}|{length}|UniProt"
#         rec = SeqRecord(
#             seq,
#             id=id,
#             description=f"{protein_name}|{ptm['ptm_note']}|{organism_name} [UniProt|Swissprot|{ptm['evidence']}]",
#         )
#         fasta_records.append(rec)
#     return fasta_records

    # with open(location, "a", encoding="UTF-8") as out:
    #     SeqIO.write(_fetch_uniprot_rest(descriptor), out, "fasta")

    SeqIO.write(_fetch_uniprot_rest(descriptor), location, "fasta")

def _fetch_uniprot_rest(descriptor):
    # url
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    url = f"{base_url}?query=reviewed:true AND ft_mod_res:{descriptor} AND existence:1&format=json"

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

            # find PTM site
            for feature in entry.get("features", []):
                if feature.get("type") != "Modified residue":
                    continue
                desc = feature.get("description", "N/A")
                if descriptor.lower() not in desc.lower():
                    continue
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                evidences = [ev.get("evidenceCode") for ev in feature.get("evidences", [])]
                eco_code = evidences[0] if evidences else "N/A"

                protein_name = f"{name}".replace(" ","__")
                organism_name = f"{organism}".replace(" ","__")
                seq = Seq(sequence)
                length = len(seq)
                id = f"{accession}|{pos}|{length}|UniProt"
                rec = SeqRecord(
                    seq,
                    id=id,
                    description=f"{protein_name}|{desc}|{organism_name} [UniProt|Swissprot|{eco_code}]",
                )
                fasta_records.append(rec)
            
        tally += len(data["results"])
        total = response.headers.get("x-total-results", "?")

        url = response.links.get("next", {}).get("url")
        if tally >= int(total):
            break

    return fasta_records
