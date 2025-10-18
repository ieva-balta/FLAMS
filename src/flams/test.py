#!/usr/bin/env python3
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import csv
import json

def query_uniprot_ptm(ptm_type, taxon_id, limit=50, output_file="uniprot_ptm_results.tsv"):
    """
    Query UniProt SPARQL endpoint for proteins with a specific PTM
    in a given organism, including gene, organism, sequence, PTM position,
    PTM description, and functional annotation.
    """

    endpoint = "https://sparql.uniprot.org/sparql"
    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)
    # sparql.setTimeout(300)

    ptm_lower = ptm_type.lower()

    # Build SPARQL query
    query = f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

    SELECT DISTINCT
    ?accession
    ?protein_name
    ?gene_name
    ?organism_name
    ?position
    ?ptm_note
    ?sequence
    ?evidence_code
    
    WHERE {{
    ?protein a up:Protein ;
            up:reviewed true ;
            up:organism taxon:{taxon_id} ;
            up:sequence ?seq .

    # Get sequence and organism name
    ?seq rdf:value ?sequence .
    ?protein up:organism ?org .
    ?org up:scientificName ?organism_name .

    # Get protein name
    OPTIONAL {{
        ?protein up:recommendedName ?recName .
        ?recName up:fullName ?protein_name .
    }}

    # Get gene name
    OPTIONAL {{
        ?protein up:encodedBy ?gene .
        ?gene skos:prefLabel ?gene_name .
    }}

    # Get PTM annotations
    ?protein up:annotation ?ptmAnn .
    ?ptmAnn a up:Modified_Residue_Annotation ;
            rdfs:comment ?ptm_note ;
            up:range/faldo:begin/faldo:position ?position .

    FILTER(CONTAINS(LCASE(?ptm_note), "{ptm_lower}"))

    # Get evidence code
    OPTIONAL {{
        ?protein up:attribution ?attr .
        ?attr up:evidence ?evidence .
        BIND(REPLACE(STR(?evidence), "^.*/(ECO_[0-9]+)$", "$1") AS ?evidence_code)
    }}

    BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession)
    }}
    ORDER BY ?accession ?position
    LIMIT {limit}
    """

    # Display query for inspection
    print("\nGenerated SPARQL Query:\n")
    print(query)
    print("\nExecuting query...\n")

    try:
        sparql.setQuery(query)
        results = sparql.query().convert()
        print(f"Query successful — {len(results['results']['bindings'])} results found.\n")
    except Exception as e:
        print(f"Query failed: {e}")
        return []

    # Parse results into structured list
    ptm_sites = []
    for row in results["results"]["bindings"]:
        site = {
            "accession": row.get("accession", {}).get("value", "N/A"),
            "protein_name": row.get("protein_name", {}).get("value", "N/A"),
            "gene_name": row.get("gene_name", {}).get("value", "N/A"),
            "organism_name": row.get("organism_name", {}).get("value", "N/A"),
            "position": row.get("position", {}).get("value", "N/A"),
            "ptm_note": row.get("ptm_note", {}).get("value", "N/A"),
            "sequence": row.get("sequence", {}).get("value", "N/A"),
            "evidence_code": row.get("evidence_code", {}).get("value", "N/A"),
            # "functionAnnotationText": row.get("functionAnnotationText", {}).get("value", "N/A")
        }
        ptm_sites.append(site)

    # Write results to TSV
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "accession",
                "protein_name",
                "gene_name",
                "organism_name",
                "position",
                "ptm_note",
                "sequence",
                "evidence_code"
                # "functionAnnotationText"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(ptm_sites)

    print(f"Results written to: {output_file}")
    return ptm_sites


if __name__ == "__main__":
    ptm = input("Enter PTM keyword (e.g. acetyl, phospho, ubiquitin): ").strip()
    taxon_input = input("Enter NCBI Taxonomy ID (e.g. 9606 for human): ").strip()
    limit_input = input("Enter number of results to retrieve (default = 50): ").strip()
    output_file = input("Enter output file name (default = uniprot_ptm_results.tsv): ").strip()

    try:
        taxon_id = int(taxon_input)
        limit = int(limit_input) if limit_input else 50
    except ValueError:
        print("Invalid input for taxon ID or limit.")
        sys.exit(1)

    if not output_file:
        output_file = "uniprot_ptm_results.tsv"

    results = query_uniprot_ptm(ptm, taxon_id, limit, output_file)

    if results:
        print("\nResults (first few shown):")
        for site in results[:5]:
            print(site)
    else:
        print("No results found.")
