#!/usr/bin/env python3
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import json
import time

def query_uniprot_ptm(ptm_type, taxon_id, limit=10, epitope=None):
    """
    Query UniProt SPARQL endpoint for proteins with a specific post-translational modification (PTM)
    in a given organism, optionally filtering by sequence epitope.
    
    Parameters:
    ptm_type (str): The type of PTM to search for (e.g., "acetyl", "phospho")
    taxon_id (int): NCBI Taxonomy ID of the organism (e.g., 9606 for human)
    limit (int): Maximum number of results to return (default 10)
    epitope (str): Optional sequence epitope to match (e.g., "VSTQ")

    Returns:
    List of dicts containing (accession, position, ptm_note)
    """
    endpoint = "https://sparql.uniprot.org/sparql"
    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)
    sparql.setTimeout(600)  # High timeout for slow PTM queries

    ptm_lower = ptm_type.lower()
    query = f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    SELECT ?accession ?position ?ptm_note
    WHERE {{
      ?protein a up:Protein .
      ?protein up:reviewed true .
      ?protein up:organism taxon:{taxon_id} .
      ?protein up:sequence ?sequence .
      ?protein up:annotation ?ann .
      ?ann a up:Modified_Residue_Annotation .
      ?ann rdfs:comment ?ptm_note .
      ?ann up:range ?range .
      ?range faldo:begin ?loc .
      ?loc faldo:position ?position .
      ?sequence rdf:value ?aaSequence .
      BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession)
      FILTER(CONTAINS(LCASE(?ptm_note), "{ptm_lower}"))
    """
    if epitope:
        query += f"""
      FILTER(SUBSTR(?aaSequence, ?position - 2, 4) = "{epitope}")
    """
    query += f"""
    }}
    ORDER BY ?accession ?position
    LIMIT {limit}
    """

    print("Generated SPARQL Query (paste into https://sparql.uniprot.org/):")
    print(query)
    print("\nExecuting...")

    try:
        start_time = time.time()
        sparql.setQuery(query)
        results = sparql.query().convert()
        elapsed_time = time.time() - start_time
        print(f"\nQuery executed in {elapsed_time:.2f} seconds")
        print("Raw results (first 500 chars for brevity):")
        print(json.dumps(results, indent=2)[:500] + "...")
        print(f"Raw bindings count: {len(results['results']['bindings'])}")
    except Exception as e:
        print(f"Query failed with exception: {e}")
        return []

    ptm_sites = []
    for row in results["results"]["bindings"]:
        ptm_sites.append({
            "accession": row["accession"]["value"],
            "position": row["position"]["value"],
            "ptm_note": row["ptm_note"]["value"]
        })

    return ptm_sites

if __name__ == "__main__":
    ptm = input("Enter PTM keyword (e.g. acetyl, phospho, ubiquitin): ").strip()
    taxon_input = input("Enter NCBI Taxonomy ID (e.g. 9606 for human): ").strip()
    limit_input = input("Enter number of results to retrieve (default = 10): ").strip()
    epitope_input = input("Enter sequence epitope (e.g. VSTQ, or press Enter for none): ").strip()
    
    try:
        taxon_id = int(taxon_input)
        limit = int(limit_input) if limit_input else 10
        epitope = epitope_input if epitope_input else None
    except ValueError:
        print("Invalid input for taxon ID or limit.")
        sys.exit(1)

    results = query_uniprot_ptm(ptm, taxon_id, limit, epitope)
    if results:
        print("\nResults:")
        for site in results:
            print(site)
    else:
        print("No results found. Steps: 1) Paste query into https://sparql.uniprot.org/. 2) Remove 'up:reviewed true' for TrEMBL. 3) Try 'phospho'. 4) Share full output.")