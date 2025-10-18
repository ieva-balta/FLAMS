#!/usr/bin/env python3
from SPARQLWrapper import SPARQLWrapper, JSON
import requests
import sys
import csv

def query_uniprot_rest_clean(ptm_type, taxon_id, limit=50):
    """
    Query UniProt REST API for proteins with PTM types based on search keyword.
    """

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    query = f'(reviewed:true) AND (organism_id:{taxon_id}) AND (cc_ptm:"*{ptm_type}*")'
    
    params = {
        'query': query,
        'format': 'json',
        'fields': 'accession,protein_name,gene_names,organism_name,sequence,ft_mod_res,cc_ptm',
        'size': limit
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=60)
        response.raise_for_status()
        
        data = response.json()
        results = data.get('results', [])
        print(f"Query successful — {len(results)} proteins found.\n")
        
        return process_rest_results_clean(results, ptm_type)
        
    except requests.exceptions.RequestException as e:
        print(f"Query failed: {e}")
        return []

def extract_ptm_details(description, ptm_keyword):
    """
    Extract additional details from PTM description.
    """
    description_lower = description.lower()
    details = []
    
    ptm_type = ptm_keyword
    
    # Extract modified amino acid
    amino_acids = ['serine', 'threonine', 'tyrosine', 'lysine', 'arginine', 'histidine', 
                   'cysteine', 'aspartate', 'glutamate', 'asparagine', 'glutamine']
    for aa in amino_acids:
        if aa in description_lower:
            details.append(aa)
            break
    
    # Extract enzymes if mentioned
    if ' by ' in description_lower:
        parts = description.split(';')
        for part in parts:
            if ' by ' in part.lower():
                details.append(part.strip())
    
    return ptm_type, '; '.join(details) if details else description

def process_rest_results_clean(results, ptm_type):
    """
    Process REST API results with PTM types based on search keyword.
    """
    ptm_sites = []
    ptm_lower = ptm_type.lower()
    
    for protein in results:
        accession = protein.get('primaryAccession', 'N/A')
        protein_name = protein.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A')
        
        # Get gene name
        gene_names = protein.get('genes', [{}])
        gene_name = 'N/A'
        if gene_names and 'geneName' in gene_names[0]:
            gene_name = gene_names[0]['geneName']['value']
        
        organism_name = protein.get('organism', {}).get('scientificName', 'N/A')
        sequence = protein.get('sequence', {}).get('value', 'N/A')
        
        # Get PTM information from modified residue features
        features = protein.get('features', [])
        for feature in features:
            if feature.get('type') == 'Modified residue':
                description = feature.get('description', '')
                if ptm_lower in description.lower():
                    position_data = feature.get('location', {}).get('start', {})
                    position = position_data.get('value', 'N/A')
                    
                    clean_ptm_type, ptm_details = extract_ptm_details(description, ptm_type)
                    
                    site = {
                        "accession": accession,
                        "protein_name": protein_name,
                        "gene_name": gene_name,
                        "organism_name": organism_name,
                        "position": position,
                        "ptm_type": clean_ptm_type,
                        "ptm_details": ptm_details,
                        "original_note": description,
                        "sequence": sequence,
                        "evidence_code": "N/A"  # From SPARQL
                    }
                    ptm_sites.append(site)
    
    # # Sort by accession and then by position (matching SPARQL's ORDER BY ?accession ?position)
    # ptm_sites.sort(key=lambda x: (
    #     x['accession'], 
    #     x['position'] if x['position'] != 'N/A' else float('inf')
    # ))
    
    return ptm_sites

def get_evidence_codes_sparql(ptm_sites, taxon_id):
    """
    Get evidence codes for all PTM sites using SPARQL.
    """
    if not ptm_sites:
        return {}
    
    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql")
    sparql.setReturnFormat(JSON)
    
    accessions = list(set(site['accession'] for site in ptm_sites))
    
    if not accessions:
        return {}
    
    accession_filter = " || ".join([f'STRENDS(STR(?protein), "/{acc}")' for acc in accessions])
    
    query = f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
    
    SELECT ?protein ?evidence ?evidence_code
    WHERE {{
      ?protein a up:Protein ;
              up:reviewed true ;
              up:organism taxon:{taxon_id} ;
              up:attribution ?att .
      
      ?att up:evidence ?evidence .
      
      BIND(REPLACE(STR(?evidence), "^.*/(ECO_[0-9]+)$", "$1") AS ?evidence_code)
      
      FILTER({accession_filter})
    }}
    """
    
    try:
        sparql.setQuery(query)
        results = sparql.query().convert()
        
        evidence_map = {}
        for row in results["results"]["bindings"]:
            protein_uri = row.get("protein", {}).get("value", "")
            evidence = row.get("evidence_code", {}).get("value", "N/A")
            
            # Extract accession from protein URI
            if "/" in protein_uri:
                accession = protein_uri.split("/")[-1]
                evidence_map[accession] = evidence
        
        print(f"Found evidence codes for {len(evidence_map)} out of {len(accessions)} proteins")
        return evidence_map
        
    except Exception as e:
        print(f"SPARQL evidence query failed: {e}")
        return {}

def query_uniprot_hybrid(ptm_type, taxon_id, limit=50, get_evidence=False):
    """
    Hybrid approach: Use REST API for speed, optionally supplement with SPARQL for evidence codes.
    """
    
    ptm_sites = query_uniprot_rest_clean(ptm_type, taxon_id, limit)
    
    if get_evidence and ptm_sites:
        evidence_map = get_evidence_codes_sparql(ptm_sites, taxon_id)
        
        # Update sites with evidence codes
        for site in ptm_sites:
            accession = site['accession']
            if accession in evidence_map:
                site['evidence_code'] = evidence_map[accession]
    
    return ptm_sites

if __name__ == "__main__":
    ptm = input("Enter PTM keyword (e.g. phospho, acetyl, ubiquitin): ").strip()
    taxon_input = input("Enter NCBI Taxonomy ID (e.g. 9606 for human): ").strip()
    limit_input = input("Enter number of results to retrieve (default = 20): ").strip()
    output_file = input("Enter output file name (default = uniprot_ptm_results.tsv): ").strip()
    evidence_choice = input("Get evidence codes? This will be slower. (y/n, default=n): ").strip().lower()

    try:
        taxon_id = int(taxon_input)
        limit = int(limit_input) if limit_input else 20
    except ValueError:
        print("Invalid input for taxon ID or limit.")
        sys.exit(1)

    if not output_file:
        output_file = "uniprot_ptm_results.tsv"
    
    get_evidence = evidence_choice == 'y'

    # Use hybrid approach
    results = query_uniprot_hybrid(ptm, taxon_id, limit, get_evidence)

    if results:
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
                    "ptm_type",
                    "ptm_details", 
                    "original_note",
                    "sequence",
                    "evidence_code"
                ],
                delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(results)

        print(f"Results written to: {output_file}")
        
        print(f"\nFirst {min(5, len(results))} results (sorted by accession → position):")
        for i, site in enumerate(results[:5], 1):
            print(f"\n{i}. {site['accession']} - {site['protein_name']}")
            print(f"   Gene: {site['gene_name']}")
            print(f"   Position: {site['position']}")
            print(f"   PTM Type: {site['ptm_type']}")
            print(f"   Evidence: {site['evidence_code']}")
            
        print(f"\nSummary: Found {len(results)} PTM sites across {len(set(s['accession'] for s in results))} proteins")
        
    else:
        print("No results found.")