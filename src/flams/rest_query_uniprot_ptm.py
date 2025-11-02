#!/usr/bin/env python3
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from fuzzywuzzy import process
import logging
import argparse
import os
import sys


# Configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("flams_ptm_fetch.log", mode="w"),
        logging.StreamHandler(sys.stdout)
    ]
)

# PTM mappings for output filtering if a PTM keyword is provided
PTM_KEYWORD_MAP = {
    "phospho": ["Phosphoserine", "Phosphothreonine", "Phosphotyrosine"],
    "acetyl": ["N6-acetyllysine"],
    "glycosyl": ["N-linked glycosylation", "O-linked glycosylation"],
    "ubiquitin": ["Ubiquitination", "Glycyl lysine isopeptide"],
    "lipid": ["S-palmitoyl cysteine"],
    "disulfide": ["Disulfide bond"],
    "sumo": ["Glycyl lysine isopeptide"]
}
VALID_PTM_KEYS = list(PTM_KEYWORD_MAP.keys())

# Valid ECO codes
VALID_ECO_CODES = {"ECO:0000269", "ECO:0000314", "ECO:0007744", "ECO:0007829"}


def validate_ptm_keyword(ptm_keyword):
    """
    Validates the user's PTM keyword and returns the list of descriptive terms 
    associated with that category for post-processing. Returns None if invalid.
    """
    ptm_keyword = ptm_keyword.strip().lower()
    for key, terms in PTM_KEYWORD_MAP.items():
        if ptm_keyword == key:
            # Return only the descriptive terms
            return [t for t in terms if len(t.split()) > 1 or t.istitle()]
    
    # Give close matches of the keyword
    matches = process.extract(ptm_keyword, VALID_PTM_KEYS, limit=3)
    logging.error("Invalid PTM keyword '%s'. Did you mean one of: %s?", ptm_keyword, ", ".join(m[0] for m in matches))
    return None


def fetch_uniprot_ptm_sequences():
    """
    Get all experimentally verified PTMs and merge features referring to the same PTM site from UniProt
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    proteins, ptm_records, seq_records = {}, {}, {}
    
    # Query using existence:1
    url = f"{base_url}?query=existence:1&format=json&size=500"

    tally = 0
    ptm_feature_types = {
        "Modified residue", "Modified residue (large scale data)", 
        "Lipidation", "Glycosylation", "Cross-link", "Disulfide bond"
    }

    while url:
        try:
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            data = r.json()
        except requests.RequestException as e:
            logging.error("UniProt API request failed: %s", str(e))
            break

        for entry in data.get("results", []):
            accession = entry.get("primaryAccession")
            sequence = entry.get("sequence", {}).get("value", "")
            protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")
            
            if not sequence: continue

            proteins[accession] = {"sequence": sequence, "protein_name": protein_name, "organism": organism}
            
            # Merge trackers
            # Dictionary to store the merged PTM entry data for this protein: 
            # Key = PTM_Description + '|' + Position
            # Value = Merged PTM entry dictionary
            merged_ptms_for_protein = {}
            # Sets to ensure evidence codes and sources are unique before merging
            
            # Process features locally
            for feature in entry.get("features", []):
                
                # Filter by PTM type
                if feature.get("type") not in ptm_feature_types:
                    continue
                
                # Get description
                desc = feature.get("description", "N/A").split(";")[0].strip()
                if not desc:
                    desc = feature.get("type", "N/A") if feature.get("type") == "Disulfide bond" else "N/A"
                
                # Get evidence and position
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                evidences = feature.get("evidences", [])
                
                # Evidence Filtering
                if not any(ev.get("evidenceCode") in VALID_ECO_CODES for ev in evidences):
                    continue
                    
                # Merge
                # Calculate unique site key -- accession is implicitly tied by iterating over entry
                site_key = f"{desc}|{pos}"
                
                # Get evidence codes and source information for the current feature
                current_eco_codes = set(ev.get("evidenceCode") for ev in evidences if ev.get("evidenceCode"))
                current_source_pairs = set(f"{ev.get('source')}:{ev.get('id')}" if ev.get('id') else ev.get('source') for ev in evidences if ev.get('source'))
                
                if site_key not in merged_ptms_for_protein:
                    # Initialize a new entry for this unique site
                    ptm_entry = {
                        "accession": accession,
                        "ptm": desc,
                        "position": pos,
                        "evidence_set": current_eco_codes,
                        "source_set": current_source_pairs,
                        "organism": organism,
                        "sequence": sequence # Store sequence temporarily for SeqRecord creation later
                    }
                    merged_ptms_for_protein[site_key] = ptm_entry
                else:
                    # Site already exists: merge evidence and source sets
                    existing_entry = merged_ptms_for_protein[site_key]
                    existing_entry["evidence_set"].update(current_eco_codes)
                    existing_entry["source_set"].update(current_source_pairs)
            
            # Finalize merged entries for this protein
            for site_key, entry in merged_ptms_for_protein.items():
                
                # Format merged sets
                eco_code_final = ";".join(sorted(entry["evidence_set"]))
                source_final = "; ".join(sorted(entry["source_set"]))

                # Final dictionary for output
                final_ptm_entry = {
                    "accession": entry["accession"],
                    "ptm": entry["ptm"],
                    "position": entry["position"],
                    "evidence": eco_code_final,
                    "source": source_final,
                    "organism": entry["organism"]
                }
                
                # Store the merged PTM record
                if entry["ptm"] not in ptm_records:
                    ptm_records[entry["ptm"]] = []
                    seq_records[entry["ptm"]] = []

                ptm_records[entry["ptm"]].append(final_ptm_entry)

                # Store the SeqRecord (using the final merged evidence/source)
                record = SeqRecord(
                    Seq(entry["sequence"]),
                    id=f"{entry['accession']}|{entry['position']}|{len(entry['sequence'])}|UniProt|{entry['ptm']}|{entry['organism']} [{eco_code_final}|{source_final}]",
                    description=""
                )
                seq_records[entry["ptm"]].append(record)

        tally += len(data.get("results", []))
        total = r.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)
        print(f"Retrieved {tally} of {total} proteins", end="\r")
        
        url = r.links.get("next", {}).get("url")
        if total != "?" and tally >= int(total): break

    return proteins, ptm_records, seq_records

def main():
    parser = argparse.ArgumentParser(description="Fetch UniProt PTM sequences (fast mode). Defaults to ALL PTMs.")
    parser.add_argument("--ptm", default=None, help="PTM keyword (e.g. phospho, acetyl). If omitted, ALL PTMs are fetched.")
    parser.add_argument("--out-dir", default="output5", help="Output directory")

    args = parser.parse_args()

    # Set PTM output mode
    if args.ptm:
        ptm_output_terms = validate_ptm_keyword(args.ptm)
        if ptm_output_terms is None:
            sys.exit(1)
        output_name = args.ptm.replace(' ', '_').lower()
        log_message = f"Fetching all PTMs, then filtering output for the '{args.ptm}' category. Duplicates merged."
    else:
        # Get all PTMs
        ptm_output_terms = None
        output_name = "all_ptms_merged"
        log_message = "Fetching and outputting ALL experimentally verified PTMs. Duplicates merged."

    os.makedirs(args.out_dir, exist_ok=True) 
    logging.info(log_message)
    print(log_message)

    try:
        # Get all PTMs with merging
        proteins, ptm_records, seq_records = fetch_uniprot_ptm_sequences()

        if not any(seq_records.values()):
            logging.warning("No PTM records found matching required evidence codes.")
            print("\n⚠ No PTM records found matching required evidence codes.")
            return

        # Post-process and filter output based on mode
        summary_file = os.path.join(args.out_dir, f"{output_name}_summary.tsv")
        total_output_entries = 0
        
        with open(summary_file, "w") as f:
            f.write("accession\tptm\tposition\tevidence\tsource\torganism\n")
            
            # Iterate through ALL PTM types fetched from UniProt
            for ptm_desc, entries in ptm_records.items():
                
                # Apply filter only if a specific PTM keyword was provided
                if ptm_output_terms is not None:
                    if not any(term.lower() in ptm_desc.lower() for term in ptm_output_terms):
                        continue
                
                # Output FASTA file for this PTM description
                safe_fasta_name = ptm_desc.replace(" ", "_").lower().replace("(", "").replace(")", "").replace("-", "_").replace("/", "_").replace(".", "")
                fasta_file = os.path.join(args.out_dir, f"{safe_fasta_name}.fasta")
                
                if seq_records[ptm_desc]:
                    # SeqRecords already contain merged evidence from fetch_uniprot_ptm_sequences
                    SeqIO.write(seq_records[ptm_desc], fasta_file, "fasta")
                    logging.info("Wrote %d sequences for %s to %s", len(seq_records[ptm_desc]), ptm_desc, fasta_file)
                
                # Output TSV file for this PTM description
                for e in entries:
                    length = len(proteins.get(e['accession'], {}).get("sequence", ""))
                    f.write(f"{e['accession']}\t{e['ptm']}\t{e['position']}\t{length}\t{e['evidence']}\t{e['source']}\t{e['organism']}\n")
                    total_output_entries += 1

        print(f"\n PTM summary written to: {summary_file}")
        logging.info("Finished. Total merged PTM entries written: %d", total_output_entries)

    except Exception as e:
        logging.error("Fatal error: %s", str(e))
        print(f"\n Error occurred: {str(e)}")

if __name__ == "__main__":
    main()