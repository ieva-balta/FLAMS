import os
import zipfile
import pandas as pd

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


# let's unzip dbptm and cplm
unzip_directory("Zenodo")
# unzip_directory("dbPTM")
# unzip_directory("CPLM")

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
                            modification = parts[4].strip()

                            df_old.append([accession, position, modification, source_name])
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

                # temp_df = pd.read_csv(file_path, sep="\t", header=None)

                # selected = temp_df.iloc[:, [1, 2, 3]].copy()

                # selected[2] = selected[2].replace("β-Hydroxybutyrylation", "b-Hydroxybutyrylation")
                # selected.columns = ["UniProt_Accession", "Position", "Modification"]

                # # Put the source 
                # selected['Source'] = source_name
                    
                # df_old.append(selected)
            
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

# Collect data from dbptm and cplm
# collect_data('dbPTM', 'dbPTM')
# collect_data('CPLM', 'CPLM')

collect_data("Zenodo")

# #Join into a single DataFrame
# if df_old:
#     combined_df = pd.concat(df_old, ignore_index=True)

#     final_df_old = combined_df.groupby(["UniProt_Accession", "Position", "Modification"])["Source"].apply(lambda x: ', '.join(sorted(set(x)))).reset_index()

#     # Save file into tsv
#     final_df_old.to_csv('PTM_combined.tsv', sep='\t', index=False)

# Create dataframe

if df_old:
    final_df_old = pd.DataFrame(df_old, columns=["UniProt_Accession", "Position", "Modification", "Source"])

    final_df_old.to_csv("PTM_old.csv", index=False)

else:
    print("No FASTA headers found")