import logging
import ssl
import pickle
import os
import pandas as pd
import requests
import time
from requests.adapters import HTTPAdapter

from urllib3 import Retry
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# checkpoint_file = "checkpoint.pkl"

# # --- Resume checkpoint if exists ---
# if os.path.exists(checkpoint_file):
#     with open(checkpoint_file, "rb") as f:
#         start_idx, records = pickle.load(f)
#     print(f"🔁 Resuming from index {start_idx}")
# else:
#     start_idx, records = 0, []

# save_counter = 0

# This is for the raw databases that they shared
# ptm_df = pd.read_csv('PTM_combined.tsv', sep='\t')
# Data from zenodo
ptm_df = pd.read_csv("PTM_old.csv")
# Load the uniprot csv
uniprot_df = pd.read_csv("Uniprot/all_records-2.0.csv")

old_accessions = set(ptm_df['UniProt_Accession'].unique())
uniprot_accessions = set(uniprot_df['Accession'].unique())

# Compute
total_unique = old_accessions.union(uniprot_accessions)

#Compute intersection
common = old_accessions.intersection(uniprot_accessions)

# PTMs unique in dbptm + CPLM
unique_old = old_accessions - uniprot_accessions

# PTMs unique in uniprot
unique_uniprot = uniprot_accessions - old_accessions

# Create a Venn diagram (requires matplotlib-venn; install with: pip install matplotlib-venn if not present)
venn2(subsets=(len(unique_uniprot), len(unique_old), len(common)),
      set_labels=('UniProt', 'CPLM/dbPTM'))
plt.title('Venn Diagram of PTMs')
plt.savefig('accessions_venn.png') #Save as PNG file
plt.close()

# Function to create a session with retries
def requests_session_with_retries():
    session = requests.Session()
    retry = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

# # Get the unique in dbptm/CPLM, the unique in uniprot and check whether those proteins are existing:1, or not
unique_old_list = list(unique_old)

session = requests_session_with_retries()
records = []
batch_size = 10000

for i in range(0, len(unique_old_list), batch_size):
    batch = unique_old_list[i:i + batch_size]
    batch_num = i // batch_size + 1
    ids = ','.join(batch)
    print(f"Processing batch {i//batch_size + 1} with {len(batch)} accessions")

    # Submit mapping job
    submit_url = "https://rest.uniprot.org/idmapping/run"
    params = {'from': 'UniProtKB_AC-ID', 'to': 'UniProtKB', 'ids': ids}
    try:
        response = session.post(submit_url, data=params)
        response.raise_for_status()
        job_id = response.json()['jobId']
        print(f"✅ Submitted job {job_id} for batch {batch_num}")
    except Exception as e:
        print(f"❌ Submission failed for batch {batch_num}: {e}")
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

            existence = result.get('to', {}).get('proteinExistence', 'Unknown')
            records.append({"Accession": acc, "Existence": existence})
            save_counter += 1
            print(f"✅ Processed {acc}: {existence}")

        for failed_acc in data.get('failedIds', []):
            records.append({"Accession": failed_acc, "Existence": "Not found"})
            save_counter += 1
            print(f"❌ {failed_acc} not found in UniProt")
        
        results_url = results_response.links.get("next", {}).get("url")
    

    if save_counter >= 100:
        with open(checkpoint_file, "wb") as f:
            pickle.dump((i + 1, records), f)
        print(f"💾 Saved checkpoint at entry {i + 1} ({len(records)} total records).")
        save_counter = 0  # reset counter
    

    time.sleep(1)

# Save final results
existence_df = pd.DataFrame(records)
existence_df.to_csv("protein_existence_summary.csv", index=False)

existence_df = pd.read_csv("protein_existence_summary.csv")

exist_1 = (existence_df["Existence"] == "1: Evidence at protein level").sum()
exist_2 = (existence_df["Existence"] == "2: Evidence at transcript level").sum()
exist_3 = (existence_df["Existence"] == "3: Inferred from homology").sum()
exist_4 = (existence_df["Existence"] == "4: Predicted").sum()
exist_5 = (existence_df["Existence"] == "5: Uncertain").sum()
exist_unknown = (existence_df["Existence"] == "Unknown").sum()
non_exist = (existence_df["Existence"] == "Not found").sum()

isoforms = existence_df["Accession"].str.contains("-").sum()

print(f"Existent proteins level 1: {exist_1}")
print(f"Existent proteins level 2: {exist_2}")
print(f"Existent proteins level 3: {exist_3}")
print(f"Existent proteins level 4: {exist_4}")
print(f"Existent proteins level 5: {exist_5}")
print(f"Existent proteins Unknown: {exist_unknown}")
print(f"Non-existent proteins: {non_exist}")
print(f"Number of isoforms: {isoforms}")