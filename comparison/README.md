# Comparison pipeline

This directory contains the scripts used to compare the UniProt datasets obtained from this FLAMS version with the datasets from dbPTM/CPLM version of [FLAMS](https://github.com/hannelorelongin/FLAMS) (database version [1.4](https://zenodo.org/records/16737546)).

The directory contains the following scripts:
- *0. dataframe_building.py* - used to prepare the two datasets for the subsequent comparisons;
- *1. df_accession_comparison.py* - used to compare the sets of accessions in both datasets;
- *2. df_comparison_unique_accessions.py* - used to investigate the set of accessions unique to the dbPTM/CPLM dataset;
- *3. df_comparison_stacked_bar.R* - used to visualize the findings from the previous step;
- *4. df_ptm_comparison.R* - used to compare the PTM record content of the datasets;
- *5. df_organism_comparison.R* - used to compare the organism diversity found in both datasets;
- *6. df_comparison_superkingdoms.py* - used to classify the species found by superkingdom;
- *7. df_comparison_updates.py* - used to investigate the effect of the next UniProt release on the dataset used in FLAMS.
