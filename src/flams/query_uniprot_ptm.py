#!/usr/bin/env python3
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import json

def query_uniprot_ptm(ptm_type, taxon_id, limit=10):
  """
  Query UniProt SPARQL endpoint for proteins with a specific post-translational modification (PTM)
  in a given organism.
  
  Parameters:
  ptm_type (str): The type of PTM to search for (e.g., "acetyl")
  taxon_id (int): NCBI Taxonomy ID of the organism (default is 9606 for human)    limit (int): Maximum number of results to return (default is 10)

  Returns:
  List of tuples containing (accession, position, ptm_label)
  """

  endpoint = "https://sparql.uniprot.org/sparql"
  sparql = SPARQLWrapper(endpoint) # UniProt SPARQL endpoint
  sparql.setReturnFormat(JSON) # return in JSON format because it's easier to parse
  sparql.setTimeout(300)

  # Query: find acetylated residues (human proteins)
  ptm_lower = ptm_type.lower()
  query = f"""
  PREFIX up: <http://purl.uniprot.org/core/>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
  PREFIX faldo: <http://biohackathon.org/resource/faldo#>
  PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  
  SELECT ?accession ?position ?ptm_note
  WHERE {{
    ?protein a up:Protein . # select proteins
    ?protein up:reviewed true . # only reviewed swiss-Prot entries
    ?protein up:organism taxon:{taxon_id} . # filter by organism
    ?protein up:sequence ?sequence . # get the sequence
    ?protein up:annotation ?ann . # get annotations
    ?ann a up:Modified_Residue_Annotation . # filter for modified residues
    ?ann rdfs:comment ?ptm_note . # get the PTM description
    ?ann up:range ?range . # get the range of the modification
    ?range faldo:begin ?loc . # get the start location
    ?loc faldo:position ?position . # get the position
    ?sequence rdf:value ?aaSequence . # get the amino acid sequence
    BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession) # extract accession from URI
    FILTER(CONTAINS(LCASE(?ptm_note), "{ptm_lower}")) # filter by PTM type
  
  }}
  ORDER BY ?accession ?position
  LIMIT {limit}
  """


  print("Generated SPARQL Query (paste into https://sparql.uniprot.org/):")
  print(query)
  print("\nExecuting...")

  # Execute the query
  try:
    # start_time = time.time()
    sparql.setQuery(query)
    results = sparql.query().convert()
    # elapsed_time = time.time() - start_time
    # print(f"\nQuery executed in {elapsed_time:.2f} seconds")
    print("Raw results (first 500 chars for brevity):")
    print(json.dumps(results, indent=2)[:500] + "...")
    print(f"Raw bindings count: {len(results['results']['bindings'])}")
  except Exception as e:
    print(f"Query failed with exception: {e}")
    return []

  # sparql.setQuery(query) # send the query to the SPARQL endpoint
  # results = sparql.query().convert() # parse the JSON results into a Python dictionary

  # Extract and return results
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
  # epitope_input = input("Enter sequence epitope (e.g. VSTQ, or press Enter for none): ").strip()
    
  try:
    taxon_id = int(taxon_input)
    limit = int(limit_input) if limit_input else 10
    # epitope = epitope_input if epitope_input else None
  except ValueError:
    print("Invalid input for taxon ID or limit.")
    sys.exit(1)

  results = query_uniprot_ptm(ptm, taxon_id, limit)
  
  if results:
    print("\nResults:")
    for site in results:
      print(site)
  else:
    print("No results found")