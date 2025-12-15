"""
@author: ieva-balta, majocava, naaattella

"""

import os
import zipfile
import pandas as pd
import urllib
import requests

"""
This script generates the dataframe files for dbPTM/CPLM database stored in Zenodo and the downloaded UniProt dataset in FLAMS. 
"""

def requests_session_with_retries():
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=3)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

# First, let's go through each directory and unzip the files

def unzip_directory(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".zip"):
            file_path = os.path.join(directory, filename)
            
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(directory)
            
            #Remove zip file
            os.remove(file_path)
            print(f"Unzipped {filename} in {directory}")


# Unzip database in Zenodo
unzip_directory("Zenodo")


df_old = []

def collect_data(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Text files
            file_path = os.path.join(root, file)

            try:
                # Read the file
                with open(file_path, 'r') as f:
                    for line in f:
                        if line.startswith(">"):
                            header = line.strip()[1:] # This is to remove >
                            parts = header.split('|')
                            accession = parts[0].strip()
                            position = parts[1].strip()
                            source_name = parts[3].strip().split(" ")[0]
                            modification = parts[4].strip().replace("__", " ")
                            organism = parts[5].strip().replace("__", " ").split("[")[0].strip()

                            df_old.append([accession, position, modification, source_name, organism])
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

# Collect data from dbptm and cplm

collect_data("Zenodo")

# Create dataframe

if df_old:
    final_df_old = pd.DataFrame(df_old, columns=["UniProt_Accession", "Position", "Modification", "Source", "Organism"])

    final_df_old.to_csv("PTM_old_organism.csv", index=False)

else:
    print("No FASTA headers found")


####################################################################
######################## UniProt Dataframe #########################
####################################################################

# Retreive data from Ebi_uniprot. Change the path if required
uniprot_path = "Ebi_uniprot/dbs_clean/dbs_clean/"

uniprot_df = []

for root, dirs, files in os.walk(uniprot_path):
    for file in files:
        if file.endswith(".fasta") or file.endswith(".fa"):
            file_path = os.path.join(root, file)

            try:
                with open(file_path, 'r') as f:
                    for line in f:
                        if line.startswith(">"):
                            header = line.strip()[1:] # Remove >
                            parts = header.split('|')
                            accession = parts[0].strip().split("_")[0]
                            position = parts[1].strip()
                            ptm_type = file.replace("-2.0.fasta", "").replace("_", " ")
                            organism = parts[5].strip().replace("__", " ").split("[")[0].strip()

                            uniprot_df.append([accession, position, ptm_type, organism])

            except Exception as e:
                print(f"Error reading {file_path}: {e}")

# Unclassified
unclassified_file = pd.read_csv("Ebi_uniprot/dbs_clean/dbs_clean/unclassified-2.0.csv")

# Function to fetch the organism
def fetch_uniprot_organism(accessions, group_size = 100):
    session = requests_session_with_retries()

    unique_accessions = list(set(accessions))
    results = {}

    group_size = 100 if len(unique_accessions)>= 100 else len(unique_accessions)

    groups = [unique_accessions[i:i + group_size] for i in range(0, len(unique_accessions), group_size)]

    print(f"Split into {len(groups)} groups for processing")
    
    for group_num, batch in enumerate(groups, 1):
        print(f"Processing group {group_num}/{len(groups)}")
        query = " OR ".join([f"accession:{acc}" for acc in batch])
        encoded_query = urllib.parse.quote(f"{query}")
        url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded_query}&format=json&size=500"

        results_response = session.get(url)
        results_response.raise_for_status()
        data = results_response.json()

        for entry in data.get("results", []):
            acc = entry["primaryAccession"]
            if acc in results:
                continue

            organism = entry["organism"]["scientificName"]
            results[acc] = organism

    return results

# Process the unclassified csv file
unc_acc = unclassified_file["Accession"].to_list()
organism_dict = fetch_uniprot_organism(unc_acc)

for _,row in unclassified_file.iterrows():
    acc = row["Accession"]
    pos = row["Position"]
    modification = row["Feature_type"].replace("_", " ")
    organism = organism_dict.get(acc)
    uniprot_df.append([acc, pos, modification, organism])

uniprot_df = pd.DataFrame(uniprot_df, columns=["UniProt_Accession", "Position", "PTM_type", "Organism"])
uniprot_df.to_csv("uniprot_organism_data_test.csv", index=False)