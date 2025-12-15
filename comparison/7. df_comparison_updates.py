"""
@author: ieva-balta, majocava, naaattella
"""

import pandas as pd
import requests

"""
This script compares the uniprot accessions with the to_remove and retain lists of accessions that will be affected in the upcoming update of UniprotKB (2026_02)
We will get the deleted/retained accessions from our UniProt dataset. 
"""


# Load the uniprot csv
uniprot_df = pd.read_csv("../uniprot_organism_data_test.csv")


# Download the files from the URLs available from UniProt's public FTP directory
url_remove = "https://ftp.ebi.ac.uk/pub/contrib/UniProt/proteomes/proteins_to_remove_from_UniProtKB.txt"
url_retain = "https://ftp.ebi.ac.uk/pub/contrib/UniProt/proteomes/proteins_retained_in_UniProtKB.txt"

def download_file(url, filename, timeout=30):
    try:
        with requests.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()

            with open(filename, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

        print(f"Downloaded {filename} successfully.")

    except requests.exceptions.RequestException as e:
        print(f"Failed to download {url}: {e}")

download_file(url_remove, "proteins_to_remove_from_UniProtKB.txt")
download_file(url_retain, "proteins_retained_in_UniProtKB.txt")

delete_list = pd.read_csv(
    "proteins_to_remove_from_UniProtKB.txt",
    sep="\t",
    header=None,
    names=["UniProt_Accession"]
)

retain_list = pd.read_csv(
    "proteins_retained_in_UniProtKB.txt",
    sep="\t",
    header=None,
    names=["UniProt_Accession"]
)


# We clean the accessions and delete the isoforms. 
cleaner_accessions = [acc.split("-")[0].strip() for acc in uniprot_df["UniProt_Accession"]]
uniprot_accessions_iso = set(cleaner_accessions)
print(f"Number of final accessions (uniprot): {len(uniprot_accessions_iso)}")

# Processing of accessions for isoforms in the to_remove dataset. 
cleaner_acccesions_d = [acc.split("-")[0].strip() for acc in delete_list["UniProt_Accession"]]
delete_list_accesion = set(cleaner_acccesions_d)
print(f"Number of accessions in the deletion list: {len(delete_list_accesion)}")

# Processing of accessions for isoforms in the retain dataset. 
cleaner_acccesions_r = [acc.split("-")[0].strip() for acc in retain_list["UniProt_Accession"]]
retain_accession = set(cleaner_acccesions_r)
print(f"Number of accessions that will remain in UniProtKb: {len(retain_accession)}")

total_unique_iso = delete_list_accesion.union(uniprot_accessions_iso)

#Compute intersection for to_remove dataset
common_iso = delete_list_accesion.intersection(uniprot_accessions_iso)

unique_uniprot_iso = uniprot_accessions_iso - common_iso

print(f"Number of accessions that will be removed from uniprot: {len(common_iso)}")
print(common_iso)

print(f"Number of accessions that will remain in uniprot: {len(unique_uniprot_iso)}")

# Compute intersection for the retained dataset. 
common_retained = retain_accession.intersection(uniprot_accessions_iso)
retain_accession_not_picked = retain_accession - uniprot_accessions_iso
accession_uniprot_removed = uniprot_accessions_iso - retain_accession

print(f"Number of accessions that will remain in UniProtKb: {len(common_retained)}")
print(f"Number of accessions that are retained in which we do not collect PTMs: {len(retain_accession_not_picked)}")
print(f"Number of accessions that that do not match the retain list: {len(accession_uniprot_removed)}")

print(accession_uniprot_removed)