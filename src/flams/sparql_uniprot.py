#!/usr/bin/env python3
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import csv
import time

class ThreadSafeSPARQL:
    def __init__(self, endpoint):
        self.endpoint = endpoint
        self._lock = threading.Lock()
    
    def query(self, query_str):
        with self._lock:
            sparql = SPARQLWrapper(self.endpoint)
            sparql.setReturnFormat(JSON)
            sparql.setMethod('POST')  # Use POST for large queries
            sparql.setQuery(query_str)
            return sparql.query().convert()

def query_uniprot_ptm_single_query(ptm_type, taxon_id, limit=50, output_file="uniprot_ptm_results.tsv"):
    """
    Single query approach: Get everything in one query - much faster
    """
    endpoint = "https://sparql.uniprot.org/sparql"
    sparql = SPARQLWrapper(endpoint)
    ptm_lower = ptm_type.lower()
    
    # Handle no limit case
    limit_clause = f"LIMIT {limit}" if limit is not None else ""
    
    print("Starting single query for PTM data with evidence codes...")
    start_time = time.time()
    
    # Single query that gets everything at once
    query = f"""
    PREFIX up: <http://purl.uniprot.org/core/>
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
            up:sequence ?seq ;
            up:organism/up:scientificName ?organism_name ;
            up:annotation ?ptmAnn .

    ?seq rdf:value ?sequence .

    OPTIONAL {{
        ?protein up:recommendedName/up:fullName ?protein_name .
    }}

    OPTIONAL {{
        ?protein up:encodedBy/skos:prefLabel ?gene_name .
    }}

    ?ptmAnn a up:Modified_Residue_Annotation ;
            rdfs:comment ?ptm_note ;
            up:range/faldo:begin/faldo:position ?position .

    # Get evidence code from protein attribution
    OPTIONAL {{
        ?protein up:attribution/up:evidence ?evidence .
        BIND(REPLACE(STR(?evidence), "^.*/(ECO_[0-9]+)$", "$1") AS ?evidence_code)
    }}

    FILTER(CONTAINS(LCASE(?ptm_note), "{ptm_lower}"))

    BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession)
    }}
    ORDER BY ?accession ?position
    {limit_clause}
    """
    
    try:
        # Use POST method for the single query
        sparql.setReturnFormat(JSON)
        sparql.setMethod('POST')
        sparql.setQuery(query)
        results = sparql.query().convert()
        
        end_time = time.time()
        print(f"Single query completed in {end_time - start_time:.2f} seconds")
        print(f"Found {len(results['results']['bindings'])} PTM sites")
        
    except Exception as e:
        print(f"Query failed: {e}")
        return []
    
    # Process results
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
            "evidence_code": row.get("evidence_code", {}).get("value", "N/A")
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
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(ptm_sites)
    
    print(f"Results written to: {output_file}")
    return ptm_sites

if __name__ == "__main__":
    ptm = input("Enter PTM keyword: ").strip()
    taxon_input = input("Enter NCBI Taxonomy ID: ").strip()
    limit_input = input("Enter number of results (default = 50, 'all' for no limit): ").strip()
    output_file = input("Enter output file name (default = uniprot_ptm_results.tsv): ").strip()

    try:
        taxon_id = int(taxon_input)
        if limit_input.lower() == 'all':
            limit = None
        else:
            limit = int(limit_input) if limit_input else 50
    except ValueError:
        print("Invalid input for taxon ID or limit.")
        sys.exit(1)

    if not output_file:
        output_file = "uniprot_ptm_results.tsv"

    # Use the single query version
    results = query_uniprot_ptm_single_query(ptm, taxon_id, limit, output_file)
    
    if results:
        unique_proteins = len(set(s['accession'] for s in results))
        print(f"\nSummary: Found {len(results)} PTM sites across {unique_proteins} proteins")
        
        # Show evidence code statistics
        evidence_counts = {}
        for site in results:
            code = site['evidence_code']
            evidence_counts[code] = evidence_counts.get(code, 0) + 1
        
        print(f"Evidence codes: {', '.join([f'{k}: {v}' for k, v in evidence_counts.items() if k != 'N/A'])}")
    else:
        print("No results found.")