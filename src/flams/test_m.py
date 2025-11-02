import re
import sys
import os
from collections import defaultdict

# --- Configuration ---
LOG_FILE_NAME = "flams_ptm_fetch.log"
OUTPUT_FILE_NAME = "log_ptm_counts.txt"

def parse_log_file():
    """Reads the log file, extracts PTM descriptions and their sequence counts."""
    ptm_counts = {}
    
    if not os.path.exists(LOG_FILE_NAME):
        print(f"Error: Log file '{LOG_FILE_NAME}' not found.")
        sys.exit(1)

    print(f"Reading log file: {LOG_FILE_NAME}")

    try:
        with open(LOG_FILE_NAME, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Check for the key logging message: "Wrote [COUNT] sequences for [PTM_NAME] to..."
                if "Wrote" in line and "sequences for" in line and "to" in line:
                    
                    # 1. Extract the count (number after "Wrote")
                    try:
                        # Split by "Wrote " to get the part containing the count and PTM name
                        count_and_ptm = line.split("Wrote ")[-1]
                        
                        # Extract the count by splitting the front part
                        count_str = count_and_ptm.split(" sequences for ")[0].strip()
                        count = int(count_str)
                        
                        # 2. Extract the PTM description (between "sequences for " and " to ")
                        ptm_desc = count_and_ptm.split(" sequences for ")[-1].split(" to ")[0].strip()

                        ptm_counts[ptm_desc] = count
                        
                    except ValueError:
                        # Skip lines where count conversion fails (shouldn't happen with standard logs)
                        continue
                    except IndexError:
                        # Skip malformed lines
                        continue

    except Exception as e:
        print(f"An error occurred during file reading: {e}")
        sys.exit(1)
        
    return ptm_counts

def write_counts_file(ptm_counts):
    """Writes the extracted PTM counts to the desired text file format."""
    
    if not ptm_counts:
        print("No PTM counts found in the log file.")
        return

    # Sort the results alphabetically by PTM description name
    sorted_counts = sorted(ptm_counts.items(), key=lambda item: item[0])

    try:
        with open(OUTPUT_FILE_NAME, 'w') as f:
            # Write header
            f.write("PTM\tCount\n")
            
            # Write data in the format PTM_DESC|COUNT
            for ptm_desc, count in sorted_counts:
                f.write(f"{ptm_desc}|{count}\n")
        
        print(f"\n PTM counts written to: {OUTPUT_FILE_NAME}")
        print(f"Total unique PTM descriptions processed: {len(ptm_counts)}")

    except Exception as e:
        print(f"Error writing output file: {e}")

if __name__ == "__main__":
    extracted_counts = parse_log_file()
    write_counts_file(extracted_counts)