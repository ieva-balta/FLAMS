#!/usr/bin/env python3
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
import csv
import threading
from concurrent.futures import ThreadPoolExecutor
import time

class ThreadSafeSPARQL:
    def __init__(self, endpoint):
        self.endpoint = endpoint
        self._lock = threading.Lock()
    
    def query(self, query_str):
        with self._lock:
            sparql = SPARQLWrapper(self.endpoint)
            sparql.setReturnFormat(JSON)
            sparql.setQuery(query_str)
            return sparql.query().convert()

def query_uniprot_ptm_parallel(ptm_type, taxon_id, limit=50, output_file="uniprot_ptm_results.tsv"):
    """
    Parallel approach: Run core query and evidence query simultaneously
    """
    endpoint = "https://sparql.uniprot.org/sparql"
    sparql = ThreadSafeSPARQL(endpoint)
    ptm_lower = ptm_type.lower()
    
    def get_core_data():
        """Get core PTM data without evidence codes"""
        print("Starting core PTM data query...")
        
        core_query = f"""
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
        
        WHERE {{
        ?protein a up:Protein ;
                up:reviewed true ;
                up:organism taxon:{taxon_id} ;
                up:sequence ?seq ;
                up:organism/up:scientificName ?organism_name .

        ?seq rdf:value ?sequence .

        OPTIONAL {{
            ?protein up:recommendedName/up:fullName ?protein_name .
        }}

        OPTIONAL {{
            ?protein up:encodedBy/skos:prefLabel ?gene_name .
        }}

        ?protein up:annotation ?ptmAnn .
        ?ptmAnn a up:Modified_Residue_Annotation ;
                rdfs:comment ?ptm_note ;
                up:range/faldo:begin/faldo:position ?position .

        FILTER(CONTAINS(LCASE(?ptm_note), "{ptm_lower}"))

        BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession)
        }}
        ORDER BY ?accession ?position
        LIMIT {limit}
        """
        
        results = sparql.query(core_query)
        print(f"Core query complete: {len(results['results']['bindings'])} results")
        return results
    
    def get_evidence_codes():
        """Get evidence codes for all proteins matching the search (runs in parallel)"""
        print("Starting evidence codes query...")
        
        # First, get all accessions that match our PTM search
        accessions_query = f"""
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        
        SELECT DISTINCT ?protein
        WHERE {{
          ?protein a up:Protein ;
                  up:reviewed true ;
                  up:organism taxon:{taxon_id} ;
                  up:annotation ?ptmAnn .
          
          ?ptmAnn a up:Modified_Residue_Annotation ;
                  rdfs:comment ?ptm_note .
          
          FILTER(CONTAINS(LCASE(?ptm_note), "{ptm_lower}"))
        }}
        LIMIT {limit}
        """
        
        accessions_result = sparql.query(accessions_query)
        proteins = [row['protein']['value'] for row in accessions_result['results']['bindings']]
        
        if not proteins:
            return {}
        
        # Get evidence codes
        protein_values = " ".join([f'<{protein}>' for protein in proteins])
        
        evidence_query = f"""
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        
        SELECT ?protein ?evidence_code
        WHERE {{
          VALUES ?protein {{ {protein_values} }}
          
          ?protein up:reviewed true ;
                   up:organism taxon:{taxon_id} ;
                   up:attribution/up:evidence ?evidence .
          
          BIND(REPLACE(STR(?evidence), "^.*/(ECO_[0-9]+)$", "$1") AS ?evidence_code)
        }}
        """
        
        evidence_results = sparql.query(evidence_query)
        
        evidence_map = {}
        for row in evidence_results["results"]["bindings"]:
            protein_uri = row.get("protein", {}).get("value", "")
            evidence_code = row.get("evidence_code", {}).get("value", "N/A")
            if "/" in protein_uri:
                accession = protein_uri.split("/")[-1]
                evidence_map[accession] = evidence_code
        
        print(f"Evidence query complete: {len(evidence_map)} codes found")
        return evidence_map
    
    # Run both queries in parallel
    print("Starting parallel queries...")
    start_time = time.time()
    
    with ThreadPoolExecutor(max_workers=2) as executor:
        # Submit both tasks
        core_future = executor.submit(get_core_data)
        evidence_future = executor.submit(get_evidence_codes)
        
        # Wait for both to complete and get results
        core_results = core_future.result()
        evidence_map = evidence_future.result()
    
    end_time = time.time()
    print(f"Parallel queries completed in {end_time - start_time:.2f} seconds")
    
    # Combine results
    ptm_sites = []
    for row in core_results["results"]["bindings"]:
        accession = row.get("accession", {}).get("value", "N/A")
        
        site = {
            "accession": accession,
            "protein_name": row.get("protein_name", {}).get("value", "N/A"),
            "gene_name": row.get("gene_name", {}).get("value", "N/A"),
            "organism_name": row.get("organism_name", {}).get("value", "N/A"),
            "position": row.get("position", {}).get("value", "N/A"),
            "ptm_note": row.get("ptm_note", {}).get("value", "N/A"),
            "sequence": row.get("sequence", {}).get("value", "N/A"),
            "evidence_code": evidence_map.get(accession, "N/A")
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

    results = query_uniprot_ptm_parallel(ptm, taxon_id, limit, output_file)
    
    if results:
        print(f"\nSummary: Found {len(results)} PTM sites across {len(set(s['accession'] for s in results))} proteins")
    else:
        print("No results found.")