import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from fuzzywuzzy import process
import logging
import argparse
import os

# ---- Logging setup ----
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ---- PTM mappings ----
PTM_KEYWORD_MAP = {
    "phospho": ["Phosphoserine", "Phosphothreonine", "Phosphotyrosine"],
    "acetyl": ["N6-acetyllysine"],
    "glycosyl": ["N-linked glycosylation", "O-linked glycosylation"],
    "ubiquitin": ["Ubiquitination"],
    "lipid": ["S-palmitoyl cysteine"],
    "disulfide": ["Disulfide bond"],
    "sumo": ["Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)"]
}
VALID_PTM_TERMS = [ptm for terms in PTM_KEYWORD_MAP.values() for ptm in terms]

# ---- Validate PTM ----
def validate_ptm_keyword(ptm_keyword):
    ptm_keyword = ptm_keyword.strip().lower()
    for key, terms in PTM_KEYWORD_MAP.items():
        if ptm_keyword == key or ptm_keyword in [t.lower() for t in terms]:
            return terms
    # Suggest close matches
    matches = process.extract(ptm_keyword, list(PTM_KEYWORD_MAP.keys()) + VALID_PTM_TERMS, limit=3)
    logging.error("Invalid PTM keyword '%s'. Did you mean one of: %s?", ptm_keyword, ", ".join(m[0] for m in matches))
    return []

def fetch_uniprot_ptm_sequences():
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/json"}
    proteins = {}
    ptm_records = {}
    seq_records = {}
    
    # Now, take all the PTM terms and fetch data, use evidence level 1 (experimental) to show that it has a protein level evidence
    # for ptm in ptm_terms:
    # url = f"{base_url}?query=reviewed:true AND existence:1&format=json&size=500" 
        # + f"AND (ft_mod_res:{ptm} OR ft_lipid:{ptm} OR ft_glyco:{ptm} OR ft_crosslnk:{ptm} OR ft_disulfide:{ptm})" 
        # 
    url = f"{base_url}?query=existence:1&format=json&size=500"

    tally = 0

    while url:
        r = requests.get(url, headers=headers)
        r.raise_for_status()
        data = r.json()

        for entry in data.get("results", []):
            accession = entry.get("primaryAccession")
            sequence = entry.get("sequence", {}).get("value", "")

            protein_name = (
                    entry.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", "")
                )
            organism = entry.get("organism", {}).get("scientificName", "")

            # Store protein info
            proteins[accession] = {
                "sequence": sequence,
                "protein_name": protein_name,
                "organism": organism
               }

            # # Initialize PTM list for this protein if not already present
            # if accession not in records:
            #     records[accession] = {
            #         "sequence": sequence,
            #         "protein_name": protein_name,
            #         "organism": organism,
            #         "ptms": []
            #         }

            #Relevant feature types:
            ptm_feature_types = {
                "Modified residue",
                "Modified residue (large scale data)",
                "Lipidation",
                "Glycosylation",
                "Cross-link",
                "Disulfide bond"
               }
            
            # PTM annotations from features
            for feature in entry.get("features", []):
                if feature.get("type") not in ptm_feature_types or feature.get("type") == "Chain":
                    continue
                    
                orig_desc = feature.get("description", "N/A")
                desc = orig_desc.split(";")[0].strip() # Clean description to get only the name
                
                if not desc or feature.get("type") == "Disulfide bond":
                    desc = feature.get("type", "N/A")

                    # if ptm.lower() not in desc.lower():
                    #     continue
                
                if desc not in ptm_records:
                    ptm_records[desc] = []
                if desc not in seq_records:
                    seq_records[desc] = []

                #Get position and evidence codes
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                evidences = feature.get("evidences", [])
                
                # Experimental evidence codes from ECO
                # ECO:0000269 (experimental evidence used in manual assertion)
                # ECO:0000314 (direct assay evidence used in manual assertion)
                # ECO:0007744 (combinatorial computational and experimental evidence used in manual assertion)
                # ECO:0007829 (combinatorial computational and experimental evidence used in automatic assertion)
                
                valid_eco_codes = {"ECO:0000269", "ECO:0000314", "ECO:0007744", "ECO:0007829", \
                                   "ECO:0000312", "ECO:0000313"}
                
                # Skip this PTM if no valid experimental evidence
                if not any(ev.get("evidenceCode") in valid_eco_codes for ev in evidences):
                    continue  

                eco_code = ";".join(ev.get("evidenceCode") for ev in evidences) if evidences else "N/A"

                # source_map = {}
                source_pairs = []
                
                for ev in evidences:
                    source_id = ev.get("source")
                    id = ev.get("id")
                    if source_id and id:
                        source_pairs.append(f"{source_id}:{id}")
                    elif source_id:
                        source_pairs.append(source_id)
                    

                        # if source not in source_map:
                        #     source_map[source] = set()
                        # if id: 
                        #     source_map[source].add(id)
                    # if ev.get("source") == "PubMed":
                    #     pubmed_ids.append(ev.get("id"))
                
                source= "; ".join(source_pairs) if source_pairs else "N/A"
                # source = ";".join([k for k in source_map.keys() if k]) if source_map else "N/A"
                # source_ids = ";".join([id for ids_set in source_map.values() for id in ids_set if id]) if source_map else "N/A"
                
                ptm_entry = {
                        "accession": accession,
                        "ptm": desc,
                        "position": pos,
                        "evidence": eco_code,
                        "source": source,
                        "organism": organism
                        # "ids": source_ids
                    }
                ptm_records[desc].append(ptm_entry)

                record = SeqRecord(
                    Seq(sequence),
                    id=f"{accession}|{pos}|{len(sequence)}|UniProt|{desc}|{organism} [{eco_code}|{source}]",
                    description = ""
                )
                
                seq_records[desc].append(record)

                #     # Aggregate PTM info for this protein
                # records[accession]["ptms"].append(f"{desc}|{pos}|{eco_code}|{source}|{source_ids}")

        tally += len(data["results"])
        total = r.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)
        print(f"Retrieved {tally} of {total} total results")

        url = r.links.get("next", {}).get("url")

        if tally >= int(total):
            break

    # Collect and display all unique PTM descriptions with counts
    print("\n List of unique PTMs found (with number of entries):")
    ptm_summary = []

    for ptm_name in sorted(ptm_records.keys()):
        count = len(ptm_records[ptm_name])
        ptm_summary.append((ptm_name, count))
        print(f"- {ptm_name}: {count} entries")

    # Save PTM names and counts to a file
    with open("ptm_list.txt", "w") as f:
        f.write("PTM\tCount\n")
        for ptm_name, count in ptm_summary:
            f.write(f"{ptm_name}\t{count}\n")

    print(f"\nPTM list written to: ptm_list.txt")

    # # Create SeqRecords with aggregated PTM info
    # seq_records = []
    # for accession, info in records.items():

    #     # Here a semicolon is better to separate multiple PTMs since it would contain similar amount of fields
    #     ptm_info = ";".join(info["ptms"])  # Combine PTM info into a single string
    #     record = SeqRecord(
    #         Seq(info["sequence"]),
    #         id=f"{accession}|UniProt",
    #         description=f"{info['protein_name']}|{ptm_info}|{info['organism']}"
    #     )
    #     seq_records.append(record)

    return proteins, ptm_records, seq_records

# ---- Main ----

# Trying to make it more general

def main():
    parser = argparse.ArgumentParser(description="Fetch UniProt protein sequences with PTMs for FLAMS database.")
    parser.add_argument("--out-dir", default="output", help="Output directory for FASTA and TSV files")
    args = parser.parse_args()
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    logging.info("Fetching sequences for all PTMs")
    print("Fetching sequences for all PTMs")
    try:
        proteins, ptm_records, seq_records = fetch_uniprot_ptm_sequences()
        if not any(seq_records.values()):
            logging.warning("No records found for any PTMs")
            print("No records found for any PTMs")
            return
        
        # Write separate files for each PTM description
        tsv_file = os.path.join(args.out_dir, "all_ptms_summary.tsv")
        with open(tsv_file, "w") as f:
            f.write("accession\tptm\tposition\tevidence\tsource\torganism\n")
        
            for ptm in sorted(ptm_records.keys()):
                
                # Clean PTM name for file naming
                ptm_safe = ptm.replace(' ', '_').lower().replace("(", "").replace(")", "").replace("-", "_").replace("/", "_").replace(".","")
                fasta_file = os.path.join(args.out_dir, f"{ptm_safe}.fasta")

                # Write FASTA file for this PTM
                if seq_records[ptm]:
                    SeqIO.write(seq_records[ptm], fasta_file, "fasta")
                    logging.info("Wrote %d sequences for %s: %s", len(seq_records[ptm]), ptm, fasta_file)
                    print(f"Wrote {len(seq_records[ptm])} sequences for {ptm}: {fasta_file}")
                else:
                    logging.warning("No sequences found for PTM %s", ptm)
                    print(f"No sequences found for PTM {ptm}")
                
                # Write TSV file for this PTM
                if ptm_records[ptm]:
                    for ptm_entry in ptm_records[ptm]:
                        length = len(proteins.get(ptm_entry['accession'], {}).get("sequence", ""))
                        f.write(f"{ptm_entry['accession']}\t{ptm_entry['ptm']}\t{ptm_entry['position']}\t{length}\t{ptm_entry['evidence']}\t{ptm_entry['source']}\t{ptm_entry['organism']}\n")
                    
                else:
                    logging.warning("No PTM records found for %s", ptm)
                    print(f"No PTM records found for {ptm}")
        
        logging.info("PTM annotation summary file written: %s", tsv_file)
        print(f"PTM annotation summary file written: {tsv_file}")
        
        # Print sample output
        print("\nSample proteins (first 3):")
        for acc, info in list(proteins.items())[:3]:
            print(f"{acc}: {info}")
        print("\nSample PTMs (first 5 per type):")
        for ptm in sorted(ptm_records.keys()):
            if ptm_records[ptm]:
                print(f"\n{ptm}:")
                for ptm_entry in ptm_records[ptm][:5]:
                    print(ptm_entry)
    except Exception as e:
        logging.error("Error processing PTMs: %s", str(e))
        print(f"Error processing PTMs: {str(e)}")

if __name__ == "__main__":
    main()