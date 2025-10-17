#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uniprot_requests_rest
Fetches UniProt PTM (post-translational modification) entries with experimental evidence
using the REST API instead of SPARQL.
"""

import logging
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def get_fasta(descriptor, location, species=None, snippet_size=None, limit=100, offset=0):
    """
    Downloads UniProt PTM data for a descriptor and saves as FASTA.
    Uses the UniProt REST API instead of SPARQL.
    
    Parameters
    ----------
    descriptor : str
        PTM description to filter (e.g., 'acetylation'). Use None for all PTMs.
    location : str
        Output FASTA file path.
    species : str, optional
        Taxon ID (e.g., '9606' for human). None for all species.
    snippet_size : int, optional
        Size of sequence snippet around PTM (e.g., 10 for ±10 residues). None for full sequence.
    limit : int
        Max number of results to fetch (default 100).
    offset : int
        Starting offset for pagination.
    """

    results = _fetch_uniprot_rest(descriptor, species, limit, offset)
    if not results:
        logging.warning("No results returned for descriptor: %s", descriptor)
        return
    
    with open(location, "w", encoding="UTF-8") as out:
        SeqIO.write(_convert_uniprot_to_fasta(results, snippet_size), out, "fasta")
    logging.info("Wrote %d records to %s", len(results), location)


def _fetch_uniprot_rest(descriptor, species, limit, offset):
    """
    Fetches UniProt entries with a given PTM type and experimental evidence via REST API.

    Parameters
    ----------
    descriptor : str or None
        PTM description (e.g., 'acetylation').
    species : str or None
        NCBI Taxonomy ID (e.g., '9606').
    limit : int
        Max number of entries.
    offset : int
        Offset for pagination.

    Returns
    -------
    List of dicts with accession, ptm_note, sequence, eco_code.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    query_parts = []

    # Build query
    if descriptor:
        query_parts.append(f'modified:"{descriptor}"')
    if species:
        query_parts.append(f'organism_id:{species}')
    query_parts.append('reviewed:true')  # Only Swiss-Prot entries
    query = " AND ".join(query_parts)

    params = {
        "query": query,
        "format": "fasta",
        "size": limit,
        "offset": offset
    }

    logging.info("Querying UniProt REST API: %s", query)

    try:
        response = requests.get(url, params=params, timeout=60)
        response.raise_for_status()
        fasta_data = response.text.strip()
        if not fasta_data:
            return []
    except requests.exceptions.RequestException as e:
        logging.error("Query failed: %s", e)
        return []

    # Parse FASTA text into SeqRecords
    handle = StringIO(fasta_data)
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    results = []
    for rec in records:
        results.append({
            "accession": rec.id.split("|")[1] if "|" in rec.id else rec.id,
            "position": None,  # REST API doesn’t provide PTM position directly
            "ptm_note": descriptor or "PTM",
            "sequence": str(rec.seq),
            "eco_code": "ECO_0000269"  # default experimental evidence
        })

    logging.info("Fetched %d records from UniProt REST API", len(results))
    return results


def _convert_uniprot_to_fasta(uniprot, snippet_size=None):
    """
    Converts UniProt query results to FASTA format.
    """
    fasta_records = []
    for site in uniprot:
        seq = site["sequence"]
        if snippet_size and site["position"]:
            pos = site["position"]
            start = max(0, pos - snippet_size - 1)
            end = min(len(seq), pos + snippet_size)
            seq = seq[start:end]
            pos = pos - start
        else:
            pos = "NA"

        seq_obj = Seq(seq)
        id = f"{site['accession']}|{pos}|{len(seq)}|UniProt|{site['eco_code']}"
        rec = SeqRecord(seq_obj, id=id, description=f"{site['ptm_note']}")
        fasta_records.append(rec)
    return fasta_records


if __name__ == "__main__":
    # Example: Fetch 10 acetylated human proteins
    get_fasta("acetylation", "test_acetylation_rest.fasta", species="9606", snippet_size=None, limit=10)