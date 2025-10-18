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
    ptms = []
    records = {}
    
    # Now, take all the PTM terms and fetch data
    # for ptm in ptm_terms:
    url = f"{base_url}?query=reviewed:true AND existence:1&format=json" 
        # + f"AND (ft_mod_res:{ptm} OR ft_lipid:{ptm} OR ft_glyco:{ptm} OR ft_crosslnk:{ptm} OR ft_disulfide:{ptm})" 
        # 
        
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

            # Initialize PTM list for this protein if not already present
            if accession not in records:
                records[accession] = {
                    "sequence": sequence,
                    "protein_name": protein_name,
                    "organism": organism,
                    "ptms": []
                    }

            #Relevant feature types:
            ptm_feature_types = {
                "Modified residue",
                "Lipidation",
                "Glycosylation",
                "Cross-link",
                "Disulfide bond"
               }
                
            # PTM annotations from features
            for feature in entry.get("features", []):
                if feature.get("type") not in ptm_feature_types:
                    continue
                    
                orig_desc = feature.get("description", "N/A")
                desc = orig_desc.split(";")[0].strip() # Clean description to get only the name
                
                if not desc:
                    desc = feature.get("type", "N/A")

                    # if ptm.lower() not in desc.lower():
                    #     continue
                pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
                evidences = feature.get("evidences", [])
                eco_code = ";".join(ev.get("evidenceCode") for ev in evidences) if evidences else "N/A"

                # Get PubMed IDS and DOIs from evidences. For this, the evidences were turned into a dictionary
                source_map = {}
                
                for ev in evidences:
                    source = ev.get("source")
                    id = ev.get("id")
                    if source:
                        if source not in source_map:
                            source_map[source] = set()
                        if id: 
                            source_map[source].add(id)
                    # if ev.get("source") == "PubMed":
                    #     pubmed_ids.append(ev.get("id"))
                    
                    #Is this necessary? Technically, the evidence do not have the DOI info, only the pubmed.
                    # if ev.get("source") == "DOI":
                    #     dois.append(ev.get("id"))

                    # citation = ref.get("citation", {})
                    # doi_id = citation.get("doi")
                    # if doi_id:
                    #     dois.append(doi_id)
                        # for ref_citation in entry.get("citation", {}).get("publication", []).get("journalCitation", []):
                        #     if ref_citation.get("idType") == "DOI":
                        #         dois.append(ref_citation.get("id", "N/A"))

                
                source = ";".join([k for k in source_map.keys() if k]) if source_map else "N/A"
                source_ids = ";".join([id for ids_set in source_map.values() for id in ids_set if id]) if source_map else "N/A"
                
                ptms.append({
                        "accession": accession,
                        "ptm": desc,
                        "position": pos,
                        "evidence": eco_code,
                        "source": source,
                        "ids": source_ids
                    })

                    # Aggregate PTM info for this protein
                records[accession]["ptms"].append(f"{desc}|{pos}|{eco_code}|{source}|{source_ids}")

        tally += len(data["results"])
        total = r.headers.get("x-total-results", "?")
        logging.info("Retrieved %d of %s total results", tally, total)
        print(f"Retrieved {tally} of {total} total results")

        url = r.links.get("next", {}).get("url")

        if tally >= int(total):
            break

    # Create SeqRecords with aggregated PTM info
    seq_records = []
    for accession, info in records.items():

        # Here a semicolon is better to separate multiple PTMs since it would contain similar amount of fields
        ptm_info = ";".join(info["ptms"])  # Combine PTM info into a single string
        record = SeqRecord(
            Seq(info["sequence"]),
            id=f"{accession}|UniProt",
            description=f"{info['protein_name']}|{ptm_info}|{info['organism']}"
        )
        seq_records.append(record)

    return proteins, ptms, seq_records


# ---- Main ----

#Trying to make it more general
def main():
    parser = argparse.ArgumentParser(description="Fetch UniProt protein sequences with all PTMs for FLAMS database.")
    parser.add_argument("--out-dir", default="output", help="Output directory for FASTA and TSV files")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Collect all PTM terms from PTM_KEYWORD_MAP
    all_ptm_terms = [ptm for terms in PTM_KEYWORD_MAP.values() for ptm in terms]
    
    logging.info("Fetching sequences for all PTMs: %s", ", ".join(all_ptm_terms))
    print(f"Fetching sequences for all PTMs: {', '.join(all_ptm_terms)}")
    
    try:
        proteins, ptms, records = fetch_uniprot_ptm_sequences()
        if not records:
            logging.warning("No records found for any PTMs")
            print("No records found for any PTMs")
            return

        # Define output file paths
        fasta_file = os.path.join(args.out_dir, "uniprot_all_ptms.fasta")
        ptm_file = os.path.join(args.out_dir, "uniprot_all_ptms.tsv")

        # Write FASTA file
        SeqIO.write(records, fasta_file, "fasta")
        logging.info("FASTA file written: %s", fasta_file)
        print(f"FASTA file written: {fasta_file}")

        # Write PTM annotations to CSV
        with open(ptm_file, "w") as f:
                f.write("accession\tptm\tposition\tevidence\tsource\tids\n")            
                for ptm_entry in ptms:
                    f.write(f"{ptm_entry['accession']}\t{ptm_entry['ptm']}\t{ptm_entry['position']}\t{ptm_entry['evidence']}\t{ptm_entry['source']}\t{ptm_entry['ids']}\n")
        logging.info("PTM annotation file written: %s", ptm_file)
        print(f"PTM annotation file written: {ptm_file}")

        # Print sample output
        print("\nSample proteins:")
        for acc, info in list(proteins.items())[:3]:
            print(f"{acc}: {info}")
        print("\nSample PTMs:")
        for ptm_entry in ptms[:5]:
            print(ptm_entry)

    except Exception as e:
        logging.error("Error processing PTMs: %s", str(e))
        print(f"Error processing PTMs: {str(e)}")

if __name__ == "__main__":
    main()

# if __name__ == "__main__":
#     ptm_input = input("Enter PTM keyword (e.g., phospho, N6-acetyllysine): ")
#     ptm_terms = validate_ptm_keyword(ptm_input)
#     if not ptm_terms:
#         exit(1)

#     # Optionally: fetch multiple PTMs in the term list
#     for ptm in ptm_terms:
#         logging.info("Fetching entries for PTM: %s", ptm)
#         proteins, ptms, records = fetch_uniprot_ptm_sequences(ptm)
#         fasta_file = f"uniprot_{ptm.replace(' ', '_').lower()}.fasta"
#         SeqIO.write(records, fasta_file, "fasta")
#         logging.info("✅ FASTA file written: %s", fasta_file)

#         # Optionally: Write PTM annotations to a separate file for FLAMS
#         ptm_file = f"uniprot_{ptm.replace(' ', '_').lower()}_ptms.csv"
#         with open(ptm_file, "w") as f:
#             f.write("accession,ptm,position,evidence\n")
#             for ptm_entry in ptms:
#                 f.write(f"{ptm_entry['accession']},{ptm_entry['ptm']},{ptm_entry['position']},{ptm_entry['evidence']}\n")
#         logging.info("✅ PTM annotation file written: %s", ptm_file)
        
#         # Example: print first 3 proteins and PTMs
#         print("\nSample proteins:")
#         for acc, info in list(proteins.items())[:3]:
#             print(acc, info)
#         print("\nSample PTMs:")
#         for ptm_entry in ptms[:5]:
#             print(ptm_entry)