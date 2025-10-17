#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uniprot_flams (PTM-specific version, FLAMS-compatible, using requests)
Fetch UniProt proteins by PTM keyword using REST API, inspired by dbPTM/CPLM.
Supports broad PTM keywords (e.g., 'phospho' for phosphorylation PTMs) and multiple species.
Produces FLAMS-compatible FASTA with sequence snippets around PTM sites.
"""
import logging
import requests
import sys
import time
import json
import os
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from urllib.parse import urlparse, parse_qs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fuzzywuzzy import process
import re
import urllib.parse

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# PTM keyword mapping
PTM_KEYWORD_MAP = {
    "phospho": ["Phosphoserine", "Phosphothreonine", "Phosphotyrosine"],
    "acetyl": ["N6-acetyllysine"],
    "glycosyl": ["N-linked glycosylation", "O-linked glycosylation"],
    "ubiquitin": ["Ubiquitination"],
    "lipid": ["S-palmitoyl cysteine"],
    "disulfide": ["Disulfide bond"],
    "sumo": ["Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)"]
}

VALID_PTM_TERMS = [
    "N6-acetyllysine", "Phosphoserine", "Phosphothreonine", "Phosphotyrosine",
    "N-linked glycosylation", "O-linked glycosylation", "Ubiquitination",
    "S-palmitoyl cysteine", "Disulfide bond",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)"
]

def validate_ptm_keyword(ptm_keyword):
    """Validate PTM keyword and suggest close matches."""
    if not ptm_keyword:
        return VALID_PTM_TERMS
    ptm_keyword = ptm_keyword.lower()
    if ptm_keyword in PTM_KEYWORD_MAP:
        return PTM_KEYWORD_MAP[ptm_keyword]

    valid_terms_lower = [t.lower() for t in VALID_PTM_TERMS]
    if ptm_keyword in valid_terms_lower:
        return [VALID_PTM_TERMS[valid_terms_lower.index(ptm_keyword)]]

    all_terms = list(PTM_KEYWORD_MAP.keys()) + VALID_PTM_TERMS
    matches = process.extract(ptm_keyword, all_terms, limit=3)
    logging.warning("Invalid PTM keyword: %s. Did you mean one of: %s?",
                    ptm_keyword, ", ".join(m[0] for m in matches))
    print(f"Error: Invalid PTM keyword '{ptm_keyword}'. Try one of: {', '.join(PTM_KEYWORD_MAP.keys())}")
    return []

def get_fasta(ptm_keyword, location, taxon_id=None, snippet_size=10, size=20):
    """Main entry point: fetch PTM proteins and save as FASTA."""
    ptm_terms = validate_ptm_keyword(ptm_keyword)
    if not ptm_terms:
        return

    results = []
    for term in ptm_terms:
        term_results = _fetch_uniprot_rest(term, taxon_id, size)
        results.extend(term_results)
        time.sleep(15)  # Delay between PTM term queries

    if not results:
        logging.warning("No results found for %s", ptm_keyword)
        print(f"No results found for '{ptm_keyword}'. Try a different PTM or taxon ID.")
        return

    with open(location, "w", encoding="UTF-8") as out:
        SeqIO.write(_convert_uniprot_to_fasta(results, snippet_size), out, "fasta")
    logging.info("✅ Wrote %d records to %s", len(results), location)
    print(f"Successfully wrote {len(results)} records to {location}")

def _fetch_uniprot_rest(ptm_term, taxon_id, size):
    """Fetch proteins with PTM annotations via UniProt REST API."""
    # Check for cached results
    cache_file = f"cache_{ptm_term}_{taxon_id or 'all'}_{size}.json"
    if os.path.exists(cache_file):
        logging.info("Loading cached results from %s", cache_file)
        with open(cache_file, "r") as f:
            return json.load(f)

    base_url = "https://rest.uniprot.org"
    query = f"reviewed:true AND annotation:{urllib.parse.quote(ptm_term)} AND existence:1"
    if taxon_id:
        query += f" AND organism_id:{taxon_id}"

    url = f"{base_url}/uniprotkb/search?query={query}&fields=accession,sequence,protein_name,ft_mod_res,organism_name&format=json&size={min(size, 25)}&sort=accession+desc"
    headers = {"Accept": "application/json", "User-Agent": "curl/8.14.1"}

    session = requests.Session()
    retries = Retry(total=3, backoff_factor=30, status_forcelist=[400, 429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))

    ptm_sites = []
    total_fetched = 0
    max_total = min(size, 500)  # Cap results
    attempt = 0
    max_attempts = 3
    delay = 10

    while url and total_fetched < max_total and attempt < max_attempts:
        logging.debug("Query URL: %s", url)
        try:
            response = session.get(url, headers=headers, timeout=120)
            if response.status_code in [400, 429, 500]:
                logging.warning("Attempt %d/%d failed with status %d for %s: %s",
                               attempt + 1, max_attempts, response.status_code, ptm_term, response.text)
                print(f"Attempt {attempt + 1}/{max_attempts} failed with status {response.status_code} for {ptm_term}: {response.text}")
                attempt += 1
                time.sleep(delay)
                delay = min(delay * 2, 60)
                continue

            response.raise_for_status()
            data = response.json()
            results = data.get("results", [])

            # Pagination handling
            url = None
            link_header = response.headers.get("Link", "")
            if link_header:
                match = re.search(r'<([^>]+)>; rel="next"', link_header)
                if match:
                    url = match.group(1)

            # Break if no results
            if not results:
                logging.info("No more results for %s. Stopping pagination.", ptm_term)
                break

            total = response.headers.get("x-total-results", "unknown")
            total_fetched += len(results)
            logging.info("REST: Got %d entries for %s (Retrieved %d of %s total results)", len(results), ptm_term, total_fetched, total)

            for entry in results:
                if total_fetched >= max_total:
                    break
                accession = entry.get("primaryAccession")
                sequence = entry.get("sequence", {}).get("value", "")
                protein_name = entry.get("proteinDescription", {}).get(
                    "recommendedName", {}
                ).get("fullName", {}).get("value", "NA").replace(" ", "__")
                species = entry.get("organism", {}).get("scientificName", "Unknown").replace(" ", "__")

                for feature in entry.get("features", []):
                    if total_fetched >= max_total:
                        break
                    if feature.get("type") != "Modified residue":
                        continue
                    description = feature.get("description", "").replace(" ", "__")
                    if ptm_term.lower() not in description.lower():
                        continue
                    evidences = feature.get("evidences", [])
                    eco_code = None
                    for ev in evidences:
                        if ev.get("evidenceCode", "") in [
                            "ECO:0000269", "ECO:0000314", "ECO:0000315",
                            "ECO:0000353", "ECO:0000316", "ECO:0000270"
                        ]:
                            eco_code = ev["evidenceCode"]
                            break
                    if not eco_code:
                        continue
                    position = feature.get("location", {}).get("start", {}).get("value", 0)
                    if position and position <= len(sequence):
                        ptm_sites.append({
                            "accession": accession,
                            "position": position,
                            "ptm_note": description,
                            "sequence": sequence,
                            "eco_code": eco_code,
                            "protein_name": protein_name,
                            "species": species
                        })
                        total_fetched += 1
            time.sleep(10)  # Delay between pagination requests
            attempt = 0  # Reset attempts after success

        except requests.exceptions.RequestException as e:
            logging.error("Attempt %d/%d failed for %s: %s", attempt + 1, max_attempts, ptm_term, e)
            print(f"Error: Attempt {attempt + 1}/{max_attempts} failed for {ptm_term}: {e}")
            attempt += 1
            if attempt >= max_attempts:
                logging.error("All retries failed for %s: %s", ptm_term, e)
                print(f"Error: UniProt query failed for {ptm_term} after {max_attempts} attempts. Try again later or use a different PTM (e.g., N6-acetyllysine).")
                return []
            time.sleep(delay)
            delay = min(delay * 2, 60)

    # Cache results
    if ptm_sites:
        with open(cache_file, "w") as f:
            json.dump(ptm_sites, f)
        logging.info("Cached results to %s", cache_file)
    return ptm_sites

def _convert_uniprot_to_fasta(uniprot, snippet_size=None):
    """Convert UniProt result dicts to FASTA, mimicking dbPTM/CPLM format."""
    fasta_records = []
    for site in uniprot:
        seq = site["sequence"]
        if not seq:
            logging.warning("Skipping entry %s: no sequence available", site["accession"])
            continue
        pos = site["position"]
        if snippet_size:
            start = max(0, pos - snippet_size - 1)
            end = min(len(seq), pos + snippet_size)
            seq = seq[start:end]
            pos = pos - start
        seq_obj = Seq(seq)
        length = len(seq)
        record_id = f"{site['accession']}|{pos}|{length}|UniProt"
        rec = SeqRecord(
            seq_obj,
            id=record_id,
            description=f"{site['protein_name']}|{site['ptm_note']}|{site['species']} [UniProt|Exp.|{site['eco_code']}]"
        )
        fasta_records.append(rec)
    return fasta_records

if __name__ == "__main__":
    ptm = input("Enter PTM keyword (e.g., phospho, N6-acetyllysine): ").strip().lower()
    taxon = input("Enter NCBI Taxonomy ID (e.g., 9606 for human, blank for all species): ").strip() or None
    limit = input("Enter number of results (1-500, default=20): ").strip()
    limit = int(limit) if limit.isdigit() else 20
    limit = min(limit, 500)  # Cap at 500
    output_file = input("Enter output FASTA file (default=uniprot_ptms.fasta): ").strip() or "uniprot_ptms.fasta"

    get_fasta(ptm, output_file, taxon_id=taxon, snippet_size=10, size=limit)