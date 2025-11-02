import requests
import pandas as pd
import logging
import sys
from typing import Dict, List, Tuple, Any

# Configuration
LOG_FILE = "ptm_fetch_df.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w"),
        logging.StreamHandler(sys.stdout)
    ]
)

# Valid ECO codes
VALID_ECO_CODES = {"ECO:0000269", "ECO:0000314", "ECO:0007744", "ECO:0007829"}

PTM_FEATURE_TYPES = {
    "Modified residue", "Modified residue (large scale data)", 
    "Lipidation", "Glycosylation", "Cross-link", "Disulfide bond"
}


def fetch_uniprot_ptm_data_and_create_df(output_pickle_file: str) -> pd.DataFrame:
    """
    Fetches all PTM data with existence:1 evidence, merges duplicates by site, 
    and returns a pandas DataFrame.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    
    # Query existence:1
    url = f"{base_url}?query=existence:1&format=json&size=500"

    tally = 0
    all_ptm_entries = []
    
    # Dictionaries to store metadata efficiently
    protein_metadata = {}

    logging.info("Starting UniProt data acquisition and merging.")
    
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
            organism = entry.get("organism", {}).get("scientificName", "")
            
            if not sequence: continue

            protein_metadata[accession] = {"sequence": sequence, "organism": organism}
            
            # Dictionary to store the merged PTM entry data for this protein: 
            # Key = PTM_Description + '|' + Position
            merged_ptms_for_protein = {}
            
            for feature in entry.get("features", []):
                
                # Filter 1: PTM Feature Type
                if feature.get("type") not in PTM_FEATURE_TYPES:
                    continue
                
                # Get Description
                desc = feature.get("description", "N/A").split(";")[0].strip()
                if not desc:
                    desc = feature.get("type", "N/A") if feature.get("type") == "Disulfide bond" else "N/A"
                
                # Get position and evidence
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                evidences = feature.get("evidences", [])
                
                # Filter ECO codes
                if not any(ev.get("evidenceCode") in VALID_ECO_CODES for ev in evidences):
                    continue
                    
                # Merge multiple features pointing to the same PTM site
                site_key = f"{desc}|{pos}"
                
                current_eco_codes = set(ev.get("evidenceCode") for ev in evidences if ev.get("evidenceCode"))
                current_source_pairs = set(f"{ev.get('source')}:{ev.get('id')}" if ev.get('id') else ev.get('source') for ev in evidences if ev.get('source'))
                
                if site_key not in merged_ptms_for_protein:
                    ptm_entry = {
                        "accession": accession, "ptm": desc, "position": pos,
                        "evidence_set": current_eco_codes, "source_set": current_source_pairs,
                        "sequence": sequence, "organism": organism
                    }
                    merged_ptms_for_protein[site_key] = ptm_entry
                else:
                    # Site already exists: merge evidence and source sets
                    existing_entry = merged_ptms_for_protein[site_key]
                    existing_entry["evidence_set"].update(current_eco_codes)
                    existing_entry["source_set"].update(current_source_pairs)
            
            # Finalize merged entries for this protein and add to master list
            for entry in merged_ptms_for_protein.values():
                
                # Format merged sets into strings
                eco_code_final = ";".join(sorted(entry["evidence_set"]))
                source_final = "; ".join(sorted(entry["source_set"]))

                all_ptm_entries.append({
                    "accession": entry["accession"],
                    "ptm_desc": entry["ptm"],
                    "position": entry["position"],
                    "sequence": entry["sequence"],
                    "evidence": eco_code_final,
                    "source": source_final,
                    "organism": entry["organism"]
                })

        tally += len(data.get("results", []))
        total = r.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)
        print(f"Retrieved {tally} of {total} proteins (Found {len(all_ptm_entries)} PTM sites)", end="\r")
        
        url = r.links.get("next", {}).get("url")
        if total != "?" and tally >= int(total): break
    
    if not all_ptm_entries:
        logging.warning("No PTM records found matching criteria.")
        return pd.DataFrame()

    df = pd.DataFrame(all_ptm_entries)
    logging.info("Data fetching complete. Total unique PTM sites found: %d", len(df))
    print(f"\nData fetching complete. Total unique PTM sites found: {len(df)}")
    
    df.to_pickle(output_pickle_file)
    logging.info("DataFrame saved to %s", output_pickle_file)
    print(f"DataFrame saved to: {output_pickle_file}")
    return df


if __name__ == "__main__":
    
    try:
        import pandas as pd
    except ImportError:
        print("Error: pandas is not installed. Please run: pip install pandas")
        sys.exit(1)

    # Output file name for the DataFrame
    OUTPUT_FILE = "uniprot_ptm_data.pkl"
    
    # Run the fetch function
    fetch_uniprot_ptm_data_and_create_df(OUTPUT_FILE)
