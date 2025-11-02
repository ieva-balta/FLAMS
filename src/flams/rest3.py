import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from fuzzywuzzy import process
import logging
import argparse
import os
import sys
import re
from typing import Dict, List, Tuple, Any


# Configuration
# Changed log filename to reflect new file
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("flams_ptm_fetch6.log", mode="w"),
        logging.StreamHandler(sys.stdout)
    ]
)


# PTM Mappings defined by REGEX patterns for efficient filtering and grouping
# Using the PTM_CATEGORY_REGEX from the last submitted working version (minus the o- ones that were empty)
PTM_CATEGORY_REGEX = {
    # Acetylation: Covers acetylation, butyrylation, propionylation, crotonylation, etc
    "acetylation": re.compile(r'(acetyl|butyryl|propionyl|crotonyl|malonyl|glutaryl|succinyl|benzoylation|lactoylation|lactylation|formylation|n-acyl|propionyl|butyryl|propionyl)', re.IGNORECASE),
    
    # ADP-Ribosylation / UMPylation / AMPylation
    "adp-ribosylation": re.compile(r'(adp-ribosyl|poly-adp-ribosyl|ribosyl|ampyl|umpyl)', re.IGNORECASE),
    
    # Amidation: Covers terminal amidation and internal modifications like phosphatidylethanolamine amidation
    "amidation": re.compile(r'(amidation|phosphatidylethanolamine amidation)', re.IGNORECASE),
    
    # Simple Acylations (Benzoylation, Butyrylation, etc.) - Grouped under Acetylation family above for broad match
    "benzoylation": re.compile(r'benzoylation', re.IGNORECASE),
    "beta-hydroxybutyrylation": re.compile(r'β-Hydroxybutyrylation|beta-hydroxybutyrylation', re.IGNORECASE),
    "butyrylation": re.compile(r'butyrylation', re.IGNORECASE),
    "crotonylation": re.compile(r'crotonylation', re.IGNORECASE),
    "lactoylation": re.compile(r'lactoylation|lactylation', re.IGNORECASE),
    "malonylation": re.compile(r'malonylation', re.IGNORECASE),
    "propionylation": re.compile(r'propionylation', re.IGNORECASE),
    "succinylation": re.compile(r'succinylation', re.IGNORECASE),

    # Biotinylation
    "biotinylation": re.compile(r'biotinylation', re.IGNORECASE),
    
    # Blocked/Modified Ends
    "blocked_amino_end": re.compile(r'blocked amino end|n-carbamoylation|pyrrolidone carboxylic acid|pyruvate', re.IGNORECASE),
    "pyrrolidone_carboxylic_acid": re.compile(r'pyrrolidone carboxylic acid', re.IGNORECASE),

    # Carboxylation/Carbamidation/Carboxyethylation/Thiocarboxylation
    "carboxylation": re.compile(r'carboxy|carbamid|carboxyethylation|gamma-carboxyglutamic acid|thiocarboxylation|n-carbamoylation|decarboxylation', re.IGNORECASE),
    
    # Lipids & Hydrophobic Attachments (Farnesyl, GPI, Myristoyl, Palmitoyl)
    "lipid": re.compile(r'(palmitoyl|myristoyl|farnesyl|geranylgeranyl|gpi-anchor|lipoylation|s-acyl|cholesterol ester|hydroxyceramide ester|s-archaeol|s-diacylglycerol|octanoylation|decanoylation|stearoylation)', re.IGNORECASE),
    "farnesylation": re.compile(r'farnesylation', re.IGNORECASE),
    "geranylgeranylation": re.compile(r'geranylgeranylation', re.IGNORECASE),
    "gpi-anchor": re.compile(r'gpi-anchor', re.IGNORECASE),
    "myristoylation": re.compile(r'myristoylation', re.IGNORECASE),
    "o-palmitoleoylation": re.compile(r'o-palmitoleoylation', re.IGNORECASE), 
    "o-palmitoylation": re.compile(r'o-palmitoylation|s-palmitoylation', re.IGNORECASE),
    "n-palmitoylation": re.compile(r'n-palmitoylation', re.IGNORECASE),
    "s-palmitoylation": re.compile(r's-palmitoylation', re.IGNORECASE),

    # Glycosylation (N-linked, O-linked, C-linked, S-linked)
    "glycosyl": re.compile(r'(glyco|glucosyl|glcnac|galnac|linked glycosylation|glycosaminoglycan|d-glucuronoylation)', re.IGNORECASE),
    "n-linked_glycosylation": re.compile(r'n-linked glycosylation', re.IGNORECASE),
    "o-linked_glycosylation": re.compile(r'o-linked glycosylation|o-linked', re.IGNORECASE),
    "c-linked_glycosylation": re.compile(r'c-linked glycosylation', re.IGNORECASE),
    "s-linked_glycosylation": re.compile(r's-linked glycosylation', re.IGNORECASE),
   
    # Hydroxylation
    "hydroxylation": re.compile(r'hydroxyl|dihydroxy', re.IGNORECASE),
    
    # Methylation
    "methylation": re.compile(r'methyl|n-methyl|di-methyl|tri-methyl|trimethyl', re.IGNORECASE),

    # Nitration
    "nitration": re.compile(r'nitr|nitro|nitrosyl|nitroso|iodination', re.IGNORECASE),

    # Oxidation/Sulfoxidation/Sulfhydration
    "oxidation": re.compile(r'oxidation|sulfoxidation|sulfhydr', re.IGNORECASE),
    
    # Phosphorylation/Dephosphorylation/Phosphoglycerylation
    "phosphorylation": re.compile(r'phospho|o-phosphate|s-phosphate|n-phosphate|dephospho|phosphoglyceryl|dietylphosphorylation', re.IGNORECASE),
    "dephosphorylation": re.compile(r'dephosphorylation', re.IGNORECASE),
    
    # Ubiquitin/SUMO/Neddylation/Pupylation
    "ubiquitination": re.compile(r'(ubiquitin|isopeptide|neddylation|pupylation|glycyl lysine|sumoylation|sumo|sentrin)', re.IGNORECASE),
    "sumoylation": re.compile(r'sumoylation', re.IGNORECASE),
    "neddylation": re.compile(r'neddylation', re.IGNORECASE),
    
    # Other specific single modifications
    "citrullination": re.compile(r'citrullination', re.IGNORECASE),
    "deamidation": re.compile(r'deamidation|deamination', re.IGNORECASE),
    "glutathionylation": re.compile(r'glutathionylation', re.IGNORECASE),
    "glycation": re.compile(r'glycation', re.IGNORECASE),
    "hmgylation": re.compile(r'hmgylation', re.IGNORECASE),
    "iodination": re.compile(r'iodination', re.IGNORECASE),
    "mgcylation": re.compile(r'mgcylation|mgylation', re.IGNORECASE),
    "serotonylation": re.compile(r'serotonylation', re.IGNORECASE),
    "sulfation": re.compile(r'sulfation', re.IGNORECASE),
    "disulfide_bond": re.compile(r'disulfide bond|cross-link|s-cyanation|s-cysteinylation|s-nitrosylation', re.IGNORECASE),
    "2-hydroxyisobutyrylation": re.compile(r'2-hydroxyisobutyrylation', re.IGNORECASE),
}

# The list of valid keys is automatically generated from the new dictionary
VALID_PTM_KEYS = list(PTM_CATEGORY_REGEX.keys())

# Valid ECO codes for experimental evidence (high quality data)
VALID_ECO_CODES = {"ECO:0000269", "ECO:0000314", "ECO:0007744", "ECO:0007829"}


def validate_ptm_keyword(ptm_keyword: str) -> bool:
    """Validates the user's PTM keyword against the PTM_CATEGORY_REGEX keys."""
    ptm_keyword = ptm_keyword.strip().lower()
    if ptm_keyword in VALID_PTM_KEYS:
        return True
    
    # Fuzzy matching for suggestions
    matches = process.extract(ptm_keyword, VALID_PTM_KEYS, limit=3)
    logging.error("Invalid PTM category '%s'. Did you mean one of: %s?", ptm_keyword, ", ".join(m[0] for m in matches))
    return False


def fetch_uniprot_ptm_data() -> Tuple[Dict, Dict, Dict]:
    """
    Fetch all PTM data with existence:1 evidence and merges duplicates by site (Accession + PTM_Description + Position)
    
    Returns (proteins dict, ptm_records dict, seq_records dict).
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    proteins, ptm_records, seq_records = {}, {}, {}
    
    # Query existence:1
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
            
            merged_ptms_for_protein = {}
            
            for feature in entry.get("features", []):
                
                # Filter PTM type
                if feature.get("type") not in ptm_feature_types:
                    continue
                
                desc = feature.get("description", "N/A").split(";")[0].strip()
                if not desc:
                    desc = feature.get("type", "N/A") if feature.get("type") == "Disulfide bond" else "N/A"
                
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                evidences = feature.get("evidences", [])
                
                # Filter experimental evidence
                if not any(ev.get("evidenceCode") in VALID_ECO_CODES for ev in evidences):
                    continue
                    
                # Merging Logic
                site_key = f"{desc}|{pos}"
                
                current_eco_codes = set(ev.get("evidenceCode") for ev in evidences if ev.get("evidenceCode"))
                current_source_pairs = set(f"{ev.get('source')}:{ev.get('id')}" if ev.get('id') else ev.get('source') for ev in evidences if ev.get('source'))
                
                if site_key not in merged_ptms_for_protein:
                    ptm_entry = {
                        "accession": accession, "ptm": desc, "position": pos,
                        "evidence_set": current_eco_codes, "source_set": current_source_pairs,
                        "organism": organism, "sequence": sequence
                    }
                    merged_ptms_for_protein[site_key] = ptm_entry
                else:
                    existing_entry = merged_ptms_for_protein[site_key]
                    existing_entry["evidence_set"].update(current_eco_codes)
                    existing_entry["source_set"].update(current_source_pairs)
            
            # Finalize and store merged entries globally
            for entry in merged_ptms_for_protein.values():
                
                # Sort sets and join with delimiters
                eco_code_final = ";".join(sorted(entry["evidence_set"]))
                source_final = "; ".join(sorted(entry["source_set"]))

                # final_ptm_entry is no longer needed for TSV output, but we keep the logic clean
                # to track the counts if necessary.
                
                ptm_desc_key = entry["ptm"] # Use the raw PTM description as the key
                
                if ptm_desc_key not in ptm_records:
                    ptm_records[ptm_desc_key] = []
                    seq_records[ptm_desc_key] = []

                # Store the simplified PTM entry (required for counting in group_and_write_ptm_files)
                ptm_records[ptm_desc_key].append(entry)

                # Store the SeqRecord using the final merged evidence/source
                record = SeqRecord(
                    Seq(entry["sequence"]),
                    id=f"{entry['accession']}|{entry['position']}|{len(entry['sequence'])}|UniProt|{ptm_desc_key}|{entry['organism']} [{eco_code_final}|{source_final}]",
                    description=""
                )
                seq_records[ptm_desc_key].append(record)

        tally += len(data.get("results", []))
        total = r.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)
        print(f"Retrieved {tally} of {total} proteins", end="\r")
        
        url = r.links.get("next", {}).get("url")
        if total != "?" and tally >= int(total): break

    return proteins, ptm_records, seq_records


def group_and_write_ptm_files(
    ptm_category: str, 
    category_regex: re.Pattern,
    proteins: Dict, 
    ptm_records: Dict, 
    seq_records: Dict,
    out_dir: str
) -> int:
    """
    Filters raw PTM data using the provided regex pattern and writes the categorized FASTA file.
    TSV output is REMOVED.
    
    Returns: The total number of PTM entries written for this category.
    """
    output_name = ptm_category.replace(' ', '_').lower()
    # summary_file line REMOVED
    fasta_file = os.path.join(out_dir, f"{output_name}.fasta") # Consolidated FASTA file
    
    total_output_entries = 0
    category_fasta_records = [] # List to accumulate ALL SeqRecords for this category
    
    logging.info("Writing files for category '%s' using regex pattern: %s", ptm_category, category_regex.pattern)

    # TSV file writing block REMOVED. Only accumulating for FASTA now.
    
    # Iterate through ALL PTM types fetched from UniProt
    for ptm_desc, entries in ptm_records.items():
        
        # Group with REGEX
        if not category_regex.search(ptm_desc):
            continue
        
        # Accumulate Sequences and Count Entries
        if ptm_desc in seq_records:
            # Accumulate FASTA records
            category_fasta_records.extend(seq_records[ptm_desc])
            # Count entries for reporting total (using the entry list size)
            total_output_entries += len(entries)
    
    # Final Consolidated FASTA
    if category_fasta_records:
        SeqIO.write(category_fasta_records, fasta_file, "fasta")
        logging.info("Wrote %d total sequences to CONSOLIDATED FASTA file: %s", len(category_fasta_records), fasta_file)
    else:
        # Added message for empty FASTA file
        print(f"\n Warning: No sequences found for category '{ptm_category}'. FASTA file was not created.")
        logging.warning("No sequences found for category '%s'. FASTA file was not created.", ptm_category)


    print(f"\n FASTA file for {ptm_category} written: {fasta_file}")
    logging.info("Finished writing files for %s. Total merged PTM entries counted: %d", ptm_category, total_output_entries)
    
    return total_output_entries


def main():
    parser = argparse.ArgumentParser(description="PTM Data Pipeline: Fetches all PTMs and outputs files grouped by category.")
    parser.add_argument("--category", default=None, help="Specific PTM category (e.g., phospho, acetyl). If omitted, ALL categories are run.")
    parser.add_argument("--out-dir", default="output_pipeline6", help="Output directory")

    args = parser.parse_args()

    # Data Acquisition
    logging.info("Starting UniProt data acquisition and merging")
    proteins, ptm_records, seq_records = fetch_uniprot_ptm_data()

    if not any(seq_records.values()):
        logging.warning("No PTM records found")
        print("\n No PTM records found")
        sys.exit(0)

    # Analysis and File Output
    
    os.makedirs(args.out_dir, exist_ok=True) 
    grand_total_entries = 0

    if args.category:
        # Targeted Run
        category = args.category.lower()
        if not validate_ptm_keyword(category):
            sys.exit(1)
        
        category_regex = PTM_CATEGORY_REGEX[category]
        logging.info("Starting targeted file output for category: %s", category)
        grand_total_entries = group_and_write_ptm_files(
            category, 
            category_regex, 
            proteins, 
            ptm_records, 
            seq_records, 
            args.out_dir
        )
    else:
        # Exhaustive Run -- loop through ALL defined categories
        logging.info("Starting exhaustive file output for all defined categories")
        
        for category, category_regex in PTM_CATEGORY_REGEX.items():
            count = group_and_write_ptm_files(
                category, 
                category_regex, 
                proteins, 
                ptm_records, 
                seq_records, 
                args.out_dir
            )
            grand_total_entries += count

    logging.info("Run Completed. Total Merged PTM Entries: %d", grand_total_entries)

if __name__ == "__main__":
    main()
