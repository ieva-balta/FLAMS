#!/usr/bin/env python3
import csv
import time
from SPARQLWrapper import SPARQLWrapper, JSON


def fetch_ptms_with_evidence(output_file="uniprot_all_ptm_filtered.tsv", batch_size=1000, max_batches=None):
    endpoint = "https://sparql.uniprot.org/sparql"
    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)

    all_results = []
    offset = 0
    batch_num = 0

    while True:
        print(f"Fetching batch {batch_num + 1} (OFFSET={offset})...")
        query = f"""
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

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
                   up:sequence ?seq ;
                   up:organism/up:scientificName ?organism_name ;
                   up:annotation ?ptmAnn .

          ?ptmAnn a up:Modified_Residue_Annotation ;
                  rdfs:comment ?ptm_note ;
                  up:range/faldo:begin/faldo:position ?position ;
                  up:evidence ?evidenceIRI .

          ?seq rdf:value ?sequence .

          OPTIONAL {{ ?protein up:recommendedName/up:fullName ?protein_name . }}
          OPTIONAL {{ ?protein up:encodedBy/skos:prefLabel ?gene_name . }}

          BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession)
          BIND(REPLACE(STR(?evidenceIRI), ".*/(ECO_\\d+)$", "$1") AS ?evidence_code)

          FILTER (
            STR(?evidenceIRI) = "http://purl.uniprot.org/evidence/evidence:1" ||
            ?evidence_code IN (
              "ECO_0000269", "ECO_0000314", "ECO_0007744",
              "ECO_0007829", "ECO_0000312", "ECO_0000313"
            )
          )
        }}
        ORDER BY ?accession ?position
        LIMIT {batch_size}
        OFFSET {offset}
        """

        sparql.setQuery(query)

        try:
            result = sparql.query().convert()
        except Exception as e:
            print(f"Error during SPARQL query: {e}")
            break

        batch = []
        for row in result["results"]["bindings"]:
            batch.append({
                "accession": row.get("accession", {}).get("value", "N/A"),
                "protein_name": row.get("protein_name", {}).get("value", "N/A"),
                "gene_name": row.get("gene_name", {}).get("value", "N/A"),
                "organism_name": row.get("organism_name", {}).get("value", "N/A"),
                "position": row.get("position", {}).get("value", "N/A"),
                "ptm_note": row.get("ptm_note", {}).get("value", "N/A"),
                "sequence": row.get("sequence", {}).get("value", "N/A"),
                "evidence_code": row.get("evidence_code", {}).get("value", "N/A")
            })

        if not batch:
            print("No more results found.")
            break

        all_results.extend(batch)
        offset += batch_size
        batch_num += 1

        if max_batches and batch_num >= max_batches:
            break

        time.sleep(1)

    if all_results:
        print(f"Writing {len(all_results)} results to {output_file}...")
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(all_results)
        print(f"Done! {len(all_results)} PTM sites across {len(set(r['accession'] for r in all_results))} proteins written.")
    else:
        print("No PTM results retrieved.")


if __name__ == "__main__":
    fetch_ptms_with_evidence("uniprot_all_ptm_filtered.tsv", batch_size=1000)
