from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import os
import sys
import argparse
import logging
import re
from typing import Dict


# Configuration
LOG_FILE = "df_to_fasta_single_file.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w"),
        logging.StreamHandler(sys.stdout)
    ]
)


# Group PTM types using regex

PTM_CATEGORY_REGEX = {
    # Acetylation & Acylations (Grouping related lysine acylations)
    "acetylation_group": re.compile(r'(acetyl|butyryl|propionyl|crotonyl|malonyl|glutaryl|succinyl|benzoylation|lactoylation|lactylation|formylation|propionyl|2-hydroxyisobutyrylation)', re.IGNORECASE),
    
    # Phosphorylation & Dephosphorylation
    "phosphorylation": re.compile(r'(phospho|o-phosphate|s-phosphate|n-phosphate|dephospho|phosphoglyceryl|dietylphosphorylation)', re.IGNORECASE),
    "dephosphorylation": re.compile(r'dephosphorylation', re.IGNORECASE),
    
    # Ubiquitin, SUMO, NEDDylation, Pupylation (UBLs and related)
    "ubl_conjugation": re.compile(r'(ubiquitin|isopeptide|neddylation|pupylation|glycyl lysine|sumoylation|sumo|sentrin)', re.IGNORECASE),
    
    # Lipids & Anchors
    "lipidation": re.compile(r'(palmitoyl|myristoyl|farnesyl|geranylgeranyl|gpi-anchor|lipoylation|s-acyl|cholesterol ester|hydroxyceramide ester|s-archaeol|s-diacylglycerol|octanoylation|decanoylation|stearoylation)', re.IGNORECASE),
    
    # Glycosylation (N-, O-, C-, S- linked and general)
    "glycosylation": re.compile(r'(glyco|glucosyl|glcnac|galnac|linked glycosylation|glycosaminoglycan|d-glucuronoylation)', re.IGNORECASE),
    
    # Structural Bonds and Cross-links
    "disulfide_bond": re.compile(r'disulfide bond|cross-link|s-cyanation|s-cysteinylation|s-nitrosylation|pyrrolylation|s-carbamoylation', re.IGNORECASE),
    
    # Simple Covalent Additions
    "amidation": re.compile(r'amidation|phosphatidylethanolamine amidation', re.IGNORECASE),
    "methylation": re.compile(r'methyl|n-methyl|di-methyl|tri-methyl|trimethyl', re.IGNORECASE),
    "hydroxylation": re.compile(r'hydroxyl|dihydroxy', re.IGNORECASE),
    "nitration": re.compile(r'nitr|nitro|nitrosyl|nitroso|iodination', re.IGNORECASE),
    "sulfation": re.compile(r'sulfation|sulfhydr|sulfoxidation|glutathionylation', re.IGNORECASE),
    "biotinylation": re.compile(r'biotinylation', re.IGNORECASE),
    "citrullination": re.compile(r'citrullination', re.IGNORECASE),
    "glycation": re.compile(r'glycation', re.IGNORECASE),
    "iodination": re.compile(r'iodination', re.IGNORECASE),
    "serotonylation": re.compile(r'serotonylation', re.IGNORECASE),
    "hmgylation": re.compile(r'hmgylation', re.IGNORECASE),
    "mgcylation_group": re.compile(r'mgcylation|mgylation', re.IGNORECASE),
    
    # Processing and Cleavage Derivatives
    "adp_ribosylation_group": re.compile(r'(adp-ribosyl|poly-adp-ribosyl|ribosyl|ampyl|umpyl)', re.IGNORECASE),
    "carboxylation_group": re.compile(r'carboxy|carbamid|carboxyethylation|gamma-carboxyglutamic acid|thiocarboxylation|n-carbamoylation|decarboxylation|pyrrolidone carboxylic acid', re.IGNORECASE),
    "deamidation_group": re.compile(r'deamidation|deamination', re.IGNORECASE),
    "oxidation": re.compile(r'oxidation', re.IGNORECASE),
    "blocked_amino_end": re.compile(r'blocked amino end|pyruvate', re.IGNORECASE),
    "pyrido_mod": re.compile(r'pyrido', re.IGNORECASE), # Catch-all for less common PTMs

    # Get everything else missed from above
    "other_modified_residue": re.compile(r'modified residue', re.IGNORECASE)
}


def get_broad_category(ptm_desc: str) -> str:
    """
    Maps a specific UniProt PTM description to a broad, single category name 
    using the defined REGEX map (First Match Wins).
    """
    # Normalize the UniProt description
    desc = ptm_desc.lower()
    
    # Iterate through the defined broad categories
    for category_name, pattern in PTM_CATEGORY_REGEX.items():
        if pattern.search(desc):
            # Return the key of the first matching category
            return category_name
    
    # Fallback for unmapped PTMs
    # This should be rare if the REGEX map is comprehensive.
    return 'unmapped_ptm_type'


def generate_fasta_files(input_pickle_file: str, output_dir: str):
    """
    Loads PTM data from a DataFrame and generates separate FASTA files, 
    grouping all PTMs of the same biological type into a single output file.
    """
    try:
        df = pd.read_pickle(input_pickle_file)
    except FileNotFoundError:
        logging.error("Input file not found: %s. Ensure fetch_and_create_df.py was run.", input_pickle_file)
        sys.exit(1)
    except Exception as e:
        logging.error("Error loading DataFrame: %s", str(e))
        sys.exit(1)

    if df.empty:
        logging.warning("DataFrame is empty. No FASTA files generated.")
        print("DataFrame is empty. No FASTA files generated.")
        return

    # Create Broad Category Column using the PTM category mapping
    df['broad_category'] = df['ptm_desc'].apply(get_broad_category)
    
    os.makedirs(output_dir, exist_ok=True)
    
    total_sites = len(df)
    unique_categories = df['broad_category'].nunique()
    logging.info("Processing %d PTM sites across %d unique broad categories.", total_sites, unique_categories)
    print(f"Processing {total_sites} PTM sites across {unique_categories} unique broad categories.")
    
    # Group by the new Broad Category
    grouped = df.groupby('broad_category')
    
    total_files_written = 0
    
    for category_name, group in grouped:
        
        # Create filename using the broad category name
        safe_name = category_name.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_").replace("/", "_").replace(".", "").lower()
        fasta_file = os.path.join(output_dir, f"{safe_name}.fasta")
        
        records = []
        
        for index, row in group.iterrows():
            
            # Create the FASTA ID line
            # Use the original 'ptm_desc' in the ID for full context
            fasta_id = (
                f"{row['accession']}|{row['position']}|{len(row['sequence'])}|UniProt|"
                f"{row['ptm_desc']}|{row['organism']} [{row['evidence']}|{row['source']}]"
            )
            
            record = SeqRecord(
                Seq(row['sequence']),
                id=fasta_id,
                description=""
            )
            records.append(record)
            
        # Write FASTA file
        if records:
            SeqIO.write(records, fasta_file, "fasta")
            logging.info("Wrote %d sequences for category '%s' to %s", len(records), category_name, fasta_file)
            total_files_written += 1
        else:
            logging.warning("Skipped writing empty file for '%s'", category_name)

    logging.info("Finished. Total FASTA files written: %d", total_files_written)
    print(f"\nFASTA generation complete. Total files written: {total_files_written}")


if __name__ == "__main__":
    
    try:
        import pandas as pd
        from Bio import SeqIO
    except ImportError:
        print("Error: pandas or Biopython is not installed. Please run: pip install pandas biopython")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Generates FASTA files from a PTM DataFrame.")
    parser.add_argument("--input-file", default="uniprot_ptm_data.pkl", help="Input pickle file from fetch_and_create_df.py")
    parser.add_argument("--out-dir", default="ptm_fasta_output_single_file", help="Output directory for FASTA files")
    
    args = parser.parse_args()
    
    generate_fasta_files(args.input_file, args.out_dir)
