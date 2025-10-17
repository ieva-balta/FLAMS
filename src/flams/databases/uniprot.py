# import requests, sys, json

# params = {
#   "query": "reviewed:true AND ft_mod_res:Phosphoserine AND existence:1",
#   "fields": [
#     "accession",
#     "protein_name",
#     "sequence",
#     "ft_mod_res"
#   ],
#   "sort": "accession desc"
# }
# headers = {
#   "accept": "application/json"
# }
# base_url = "https://rest.uniprot.org/uniprotkb/search"

# response = requests.get(base_url, headers=headers, params=params)

# tally = 0
# while response:
#     r=get_url(response)
#     total = r.headers.get("x-total-results")
#     print(f"Retrieved {tally} of {total} total results")
#     tally += len(r.get("results", []))
#     url = r.links.get("next", {}).get("url")

# if not response.ok:
#   response.raise_for_status()
#   sys.exit()

# data = response.json()

# # Extract what we care about
# for entry in data.get("results", []):
#     accession = entry.get("primaryAccession")
#     seq_info = entry.get("sequence", {})
#     sequence = seq_info.get("value", "")
#     protein_name = entry.get("proteinDescription", {}).get(
#         "recommendedName", {}
#     ).get("fullName", {}).get("value", "N/A")

#     print(f"\nAccession: {accession}")
#     print(f"Protein name: {protein_name}")
#     print(f"Sequence length: {len(sequence)}")

#     for feature in entry.get("features", []):
#         if feature.get("type") != "Modified residue":
#             continue
#         desc = feature.get("description", "N/A")
#         pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
#         evidences = [ev.get("evidenceCode") for ev in feature.get("evidences", [])]
#         print(f"  - {desc} at position {pos} (Evidence: {', '.join(evidences) or 'None'})")

# # print(json.dumps(data, indent=2))



# # Code to paginate through all results according to the tutorial
# url = f"{base_url}/uniprotkb/search?query=reviewed:true AND ft_mod_res:Phosphoserine AND existence:1"

# tally = 0
# while url:
#     r = requests.get(url, headers=headers)
#     data = r.json()
#     total = r.headers.get("x-total-results")
#     print(f"Retrieved {tally} of {total} total results")
#     tally += len(data["results"])
#     url = r.links.get("next", {}).get("url")
import requests
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from fuzzywuzzy import process
import logging

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

def fetch_uniprot_ptm_sequences(ptm_term):
  base_url = "https://rest.uniprot.org/uniprotkb/search"
  headers = {"Accept": "application/json"}

  # Query for reviewed entries with specified PTM and protein existence level 1
  url = f"{base_url}?query=reviewed:true AND ft_mod_res:{ptm_term} AND existence:1&format=json"

  proteins = {}
  ptms = []
  records = []
  tally = 0

  while url:
      r = requests.get(url, headers=headers)
      r.raise_for_status()
      data = r.json()

      # Extract entries
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
          
          proteins[accession] = {
              "sequence": sequence,
              "protein_name": protein_name,
              "organism": organism
          }

          #PTM annotations
          for feature in entry.get("features", []):
              if feature.get("type") != "Modified residue":
                  continue
              desc = feature.get("description", "N/A")
              if ptm_term.lower() not in desc.lower():
                  continue
              pos = feature.get("location", {}).get("start", {}).get("value", "N/A")
              evidences = [ev.get("evidenceCode") for ev in feature.get("evidences", [])]
              eco_code = evidences[0] if evidences else "N/A"

              ptms.append({
                  "accession": accession,
                  "ptm": desc,
                  "position": pos,
                  "evidence": eco_code
              })
          
              record = SeqRecord(
                Seq(sequence),
                id=f"{accession}|{pos}|UniProt",
                description=f"{protein_name}|{desc}|{pos}|{eco_code}|{organism}"
              )
              records.append(record)

      tally += len(data["results"])
      total = r.headers.get("x-total-results", "?")
      logging.info("Retrieved %d of %s total results", tally, total)
      print(f"Retrieved {tally} of {total} total results")

      # Follow pagination link
      url = r.links.get("next", {}).get("url")

      if tally >= int(total):
          break
  
  return proteins, ptms, records

# print(f"Total records fetched: {len(records)}")
# SeqIO.write(records, "uniprot_phosphoserine.fasta", "fasta")
# print("✅ FASTA file written: uniprot_phosphoserine.fasta")


# ---- Main ----
if __name__ == "__main__":
    ptm_input = input("Enter PTM keyword (e.g., phospho, N6-acetyllysine): ")
    ptm_terms = validate_ptm_keyword(ptm_input)
    if not ptm_terms:
        exit(1)

    # Optionally: fetch multiple PTMs in the term list
    for ptm in ptm_terms:
        logging.info("Fetching entries for PTM: %s", ptm)
        proteins, ptms, records = fetch_uniprot_ptm_sequences(ptm)
        fasta_file = f"uniprot_{ptm.replace(' ', '_').lower()}.fasta"
        SeqIO.write(records, fasta_file, "fasta")
        logging.info("✅ FASTA file written: %s", fasta_file)

        # Example: print first 3 proteins and PTMs
        print("\nSample proteins:")
        for acc, info in list(proteins.items())[:3]:
            print(acc, info)
        print("\nSample PTMs:")
        for ptm_entry in ptms[:5]:
            print(ptm_entry)