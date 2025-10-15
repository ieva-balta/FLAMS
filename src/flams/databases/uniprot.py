#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" uniprot
functions to fetch and parse uniprot entries
based on scripts for cplm4 and dbptm
"""

from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_fasta(descriptor, location):
    """
    will download uniprot data for each descriptor (modification) seperately
    and save it in a fasta format
    return the multi fasta file for the specified decriptor
    """
    
    with open(location, "a", encoding="UTF-8") as out:
        SeqIO.write(_convert_uniprot_to_fasta(_fetch_uniprot(descriptor)), out, "fasta")


### need to figure out what to do with different PTM names
def _fetch_uniprot(descriptor):
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

def _convert_uniprot_to_fasta(uniprot):
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