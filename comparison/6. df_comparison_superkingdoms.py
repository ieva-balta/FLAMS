
"""
@author: ieva-balta, majocava, naaattella
"""

import pickle
import os
import pandas as pd
import requests
import time
from requests.adapters import HTTPAdapter
import urllib.parse
from Bio import Entrez
import time
from urllib3 import Retry
import matplotlib.pyplot as plt

"""
This script will retrieve the superkingdom entries from NCBI and/or Uniprot REST API.
This is to compare the organism data. 
Old refers to the dbPTM/CPLM database.
"""

# Set your email for NCBI Entrez (required)
Entrez.email = "majocv97@gmail.com" # Change to email

# Function to get superkingdom from species name using NCBI taxonomy (This is for dbPTM/CPLM)
def get_superkingdom(species):
    try:
        # Format species, species lowercase
        species_parts = species.split()
        if len(species_parts)>=2:
            species = f"{species_parts[0].capitalize()} {species_parts[1].lower()}"
        else:
            species = species.capitalize()
        
        # Search for taxonomy ID
        handle = Entrez.esearch(db="taxonomy", term=species)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            return 'Unknown'
        tax_id = record["IdList"][0]

        # Fetch taxonomy details
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        data = Entrez.read(handle)
        handle.close()

        if not data:
            return "Unknown"
        
        # Get lineage
        lineage = data[0].get('Lineage', '')
        if not lineage:
            return 'Unknown'
        
        levels = [level.strip() for level in lineage.split(";")]

        # Determine superkingdowm
        if levels and levels[0] == 'cellular organisms' and len(levels) >1:
            superkingdom = levels[1]
        elif "Viruses" in levels:
            superkingdom = "Viruses"
        else:
            superkingdom = levels[0] if levels else 'Unknown'
        
        return superkingdom.capitalize()

    except Exception as e:
        print(f"Error for {species}: {e}")
        return "Unknown"
    
def get_superkinddom_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    r = requests.get(url)

    data = r.json()
    lineage = data["organism"].get("lineage", [])
    if not lineage:
        return "Unknown"
    
    return lineage[0].capitalize()

# Data from zenodo
ptm_old = pd.read_csv("../PTM_old_organism.csv")
# Load the uniprot csv
uniprot_df = pd.read_csv("../uniprot_organism_data_test.csv")

# Fill specific missing Organism values based on known species
uniprot_df.loc[uniprot_df['UniProt_Accession'] == 'A0A5S8WF70', 'Organism'] = 'Leptolyngbya sp. JSC-1'
uniprot_df.loc[uniprot_df['UniProt_Accession'] == 'Q845W8', 'Organism'] = 'Pseudomonas dacunhae'


# Funtion to clean the species names. 
def clean_species(x):
    if pd.isna(x) or (isinstance(x, str) and x.strip() == ""):
        return pd.NA
    
    if not isinstance(x, str):
        return pd.NA
    
    x = x.lower().strip()
    x = x.split('(')[0].strip()
    
    parts = x.split()
    if not parts:
        return pd.NA
    
    if parts[0].startswith("recombinant"):
        parts = parts[1:]
    
    parts = parts[:min(2, len(parts))]
    
    return ' '.join(parts)

# Apply to both DataFrames
uniprot_df['Species'] = uniprot_df['Organism'].apply(clean_species)
ptm_old['Species'] = ptm_old['Organism'].apply(clean_species)

# Find if any accession don't have organism data. 
bad = uniprot_df[ptm_old["Organism"].isna()]
print(f"Rows with missing organisms:\n", bad)

# Compare the organism distributions
old_acc = set(ptm_old['UniProt_Accession'])
uniprot_acc = set(uniprot_df['UniProt_Accession'])

shared_acc = old_acc.intersection(uniprot_acc)
unique_old_acc = old_acc - uniprot_acc
unique_uniprot_acc = uniprot_acc - old_acc

# Which organisms are present both datasets, unique to old, unique to uniprot
organism_shared = set(ptm_old[ptm_old['UniProt_Accession'].isin(shared_acc)]['Species']).union(
    set(uniprot_df[uniprot_df['UniProt_Accession'].isin(shared_acc)]['Species'])
)
organism_old_only = set(ptm_old[ptm_old['UniProt_Accession'].isin(unique_old_acc)]['Species'])
organism_uniprot_only = set(uniprot_df[uniprot_df['UniProt_Accession'].isin(unique_uniprot_acc)]['Species'])

species_old = set(ptm_old['Species'])
species_uniprot = set(uniprot_df['Species'])

organism_shared = species_old.intersection(species_uniprot)
organism_old_only = species_old - species_uniprot
organism_uniprot_only = species_uniprot - species_old

print(f"Organisms in both datasets: {len(organism_shared)}")

print(f"\nOrganisms unique to old dataset: {len(organism_old_only)}")

print(f"\nOrganisms unique to UniProt dataset: {len(organism_uniprot_only)}")

# Get unique species
unique_species_uniprot = uniprot_df["Species"].unique()

# species to accesions
species_to_accessions = {}
unique_accessions = set()
for sp in unique_species_uniprot:
    accessions = uniprot_df[uniprot_df["Species"] == sp]["UniProt_Accession"].unique().tolist()
    species_to_accessions[sp] = accessions
    unique_accessions.update(accessions)
    

# Function to create a session with retries
def requests_session_with_retries():
    session = requests.Session()
    retry = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

# Get the unique UniPrOT ORGANISM
unique_uniprot_list = list(uniprot_acc)

session = requests_session_with_retries()
records = []
batch_size = 10000
save_counter = 0
checkpoint_file = "checkpoint.pkl"

for i in range(0, len(unique_uniprot_list), batch_size):
    batch = unique_uniprot_list[i:i + batch_size]
    batch_num = i // batch_size + 1
    ids = ','.join(batch)
    print(f"Processing batch {i//batch_size + 1} with {len(batch)} accessions out of {batch_num}")

    # Submit mapping job
    submit_url = "https://rest.uniprot.org/idmapping/run"
    params = {'from': 'UniProtKB_AC-ID', 'to': 'UniProtKB', 'ids': ids}
    try:
        response = session.post(submit_url, data=params)
        response.raise_for_status()
        job_id = response.json()['jobId']
        print(f"Submitted job {job_id} for batch {batch_num}")
    except Exception as e:
        print(f"Submission failed for batch {batch_num}: {e}")
        continue

    # Poll job until finished
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status_response = session.get(status_url)
        status_response.raise_for_status()
        status_json = status_response.json()
        print(f"Checking status for job {job_id}: {status_json.get('jobStatus')}.")
        if status_json.get('jobStatus') is None:
            print(f"Job {job_id} finished. Fetching results")
            break
        if status_json.get('jobStatus') in ['FAILED', 'ERROR']:
            print(f"Job failed: {job_id}")
            break

        time.sleep(5)

    # Fetch results
    results_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}?size=500"
    
    
    while results_url:
        results_response = session.get(results_url)
        results_response.raise_for_status()
        data = results_response.json()
        for result in data.get('results', []):
            acc = result['from']
            to_data = result.get('to')

            primary_acc = to_data.get("primaryAccession", acc)

            scientific_name = to_data.get("organism", {}).get("scientificName", "Unknown")
            lineage = to_data.get("organism", {}).get("lineage", [])
            superkingdom = lineage[0].capitalize() if lineage else "Unknown"
            subdivision = lineage[1].capitalize() if len(lineage) > 1 else "Unknown"
            
        
            print(f"Processed {acc}: {superkingdom}, {subdivision}")

            records.append({
                "accession" : acc,
                "primary_accession": primary_acc,
                "scientific_name": scientific_name,
                "superkingdom": superkingdom,
                "subdivision" : subdivision})

            save_counter += 1

        results_url = results_response.links.get("next", {}).get("url")
    

    if save_counter >= 100:
        with open(checkpoint_file, "wb") as f:
            pickle.dump((i + 1, records), f)
        print(f"💾 Saved checkpoint at entry {i + 1} ({len(records)} total records).")
        save_counter = 0  # reset counter
    

    time.sleep(1)

# Save final results
existence_df = pd.DataFrame(records)
existence_df.to_csv("superkingdom_unique_uniprot.csv", index=False)

##################################################
########## Unique species from ptm_old ###########
##################################################

species_to_sk_old = {}
count = 0
for sp in organism_old_only:
    count += 1
    print(f"Processing {count} species unique to old database")
    superkingdom = get_superkingdom(sp)
    
    # Get the accession
    if superkingdom == "Unknown":
        accessions = ptm_old[ptm_old["Species"] == sp]["UniProt_Accession"].unique()

        if len(accessions) > 0:
            accession = accessions[0]
            superkingdom = get_superkinddom_uniprot(accession)

    species_to_sk_old[sp] = superkingdom
    time.sleep(0.4) # delay for rate limits


ptm_old["superkingdom"] = ptm_old["Species"].map(species_to_sk_old)

# save updated ptm_old
ptm_old.to_csv("ptm_old_with_superkingdom_unique.csv", index = False)

# Make comparison 
existence_old_df = pd.read_csv("ptm_old_with_superkingdom_unique.csv")
existence_df = pd.read_csv("superkingdom_unique_uniprot.csv")

unique_uniprot_species_df = uniprot_df[uniprot_df['Species'].isin(organism_uniprot_only)]

# First, merge superkingdom back to uniprot_df
unique_uniprot_species_df = unique_uniprot_species_df.merge(
    existence_df[['accession', 'superkingdom']],
    left_on='UniProt_Accession',
    right_on='accession',
    how='left'
).drop(columns='accession')

# Get unique species with their superkingdom
species_superkingdom_unique = unique_uniprot_species_df.groupby('Species')['superkingdom'].first().reset_index()

# Count number of unique species per superkingdom
species_counts_unique = species_superkingdom_unique.groupby('superkingdom')['Species'].nunique().reset_index()

# Rename columns for clarity
species_counts_unique.columns = ['superkingdom', 'unique_species_count']

# Sort by count
species_counts_unique = species_counts_unique.sort_values('unique_species_count', ascending=False)

print(species_counts_unique)

# Verify totals
total_unique_species = len(organism_uniprot_only)
print(f"Total unique species in UniProt-unique organisms: {total_unique_species}")
total_counted = species_counts_unique['unique_species_count'].sum()
print(f"Total counted in superkingdoms: {total_counted}")

# Filter rows for 'Unclassified sequences' from the DataFrame
unclassified_df = unique_uniprot_species_df[unique_uniprot_species_df['superkingdom'] == 'Unclassified sequences']

# Group by 'Species' and collect unique accessions as lists
unclassified_species_accessions = unclassified_df.groupby('Species')['UniProt_Accession'].unique().apply(list).reset_index()

# Print the result as a table
print("Unclassified sequences species and their accessions:")
print(unclassified_species_accessions)

########################################
# For the case of the shared organisms #
########################################

shared_species_df = uniprot_df[uniprot_df['Species'].isin(organism_shared)]

# First, merge superkingdom back to uniprot_df
shared_species_df = shared_species_df.merge(
    existence_df[['accession', 'superkingdom']],
    left_on='UniProt_Accession',
    right_on='accession',
    how='left'
).drop(columns='accession')

# Get unique species with their superkingdom
species_superkingdom_shared = shared_species_df.groupby('Species')['superkingdom'].first().reset_index()

# Count number of unique species per superkingdom
species_counts_shared = species_superkingdom_shared.groupby('superkingdom')['Species'].nunique().reset_index()

# Rename columns for clarity
species_counts_shared.columns = ['superkingdom', 'shared_species_count']

# Sort by count
species_counts_shared = species_counts_shared.sort_values('shared_species_count', ascending=False)

print(species_counts_shared)

# Filter and print the species in 'Unclassified sequences'
unclassified_species_shared = species_superkingdom_shared[species_superkingdom_shared['superkingdom'] == 'Unclassified sequences']['Species'].tolist()
print("Unclassified sequences species:")
print(unclassified_species_shared)

############################################
############ PTM Old species ###############
############################################

old_species_df = existence_old_df[existence_old_df['Species'].isin(organism_old_only)]

# Similar counting for old unique
species_superkingdom_old = old_species_df.groupby('Species')['superkingdom'].first().reset_index()
species_counts_old = species_superkingdom_old.groupby('superkingdom')['Species'].nunique().reset_index()
species_counts_old.columns = ['superkingdom', 'unique_species_count']
species_counts_old = species_counts_old.sort_values('unique_species_count', ascending=False)
print("Counts for unique old species:")
print(species_counts_old)

# Unclassified for old
unclassified_old = old_species_df[old_species_df['superkingdom'] == 'Unknown']
unclassified_species_accessions_old = unclassified_old.groupby("Species")["UniProt_Accession"].unique().apply(list).reset_index()

# Print the result as a table
print("Unclassified sequences species and their accessions in old unique:")
print(unclassified_species_accessions_old)