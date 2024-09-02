import numpy as np
import pandas as pd
import os
import argparse
from config_and_print import resolutions, chromosomes, chrom_file

# Print initial message
print("Starting process_chromsizes.py script...")

# Ensure resolutions is treated as a tuple or list of strings
if isinstance(resolutions, str):
    resolutions = (resolutions,)

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]

def parse_resolution(resolution_str):
    if ':' in resolution_str:
        resolution_value, resolution_label = resolution_str.split(':')
        try:
            resolution = int(resolution_value)
            return resolution, resolution_label
        except ValueError:
            raise ValueError(f"Resolution value should be an integer: '{resolution_value}' in '{resolution_str}'")
    else:
        raise ValueError(f"Invalid resolution format: '{resolution_str}'. Expected format 'value:label', e.g., '1000000:1Mb'.")

resolution, resolution_label = parse_resolution(resolution_str)

# Load chromosome sizes from a file specified by the chrom_file variable
try:
    print(f"Loading chromosome sizes from: {chrom_file}")
    chromsize = pd.read_csv(chrom_file, sep='\t', header=None, index_col=0).to_dict()[1]
    print(f"Loaded chromosome sizes: {chromsize}")
except FileNotFoundError:
    print(f"File not found: {chrom_file}")
    raise
except Exception as e:
    print(f"Error loading chromosome sizes: {e}")
    raise

# Set up argument parsing for just the directory
parser = argparse.ArgumentParser(description="Process chromosome sizes and create BED files.")
parser.add_argument('dir', help='Existing output directory for BED files')

args = parser.parse_args()

# Print the output directory
print(f"Output directory: {args.dir}")

# Loop over each resolution and chromosome to create BED segments
for res_str, label in resolutions.items():
    print(f"Processing resolution: {res_str}, label: {label}")
    res = int(res_str)  # Convert resolution string to integer
    for c in chromsize:
        print(f"Processing chromosome: {c}, resolution: {res}")
        # Calculate the number of segments based on chromosome size and resolution
        ngene = int(chromsize[c] // res) + 1
        print(f"Number of segments (ngene): {ngene}")
        # Generate BED format data for each segment
        bed = [[c, i*res, (i+1)*res] for i in range(ngene)]
        # Ensure the last segment extends to the end of the chromosome
        bed[-1][-1] = chromsize[c]
        print(f"Generated BED data for chromosome {c}: {bed[-1]}")

        # Save the BED data to a file, one for each chromosome
        outdir = f'{args.dir}/bins'
        print(f"Checking if directory exists: {outdir}")
        if not os.path.exists(outdir):
            print(f"Directory {outdir} does not exist. Raising FileNotFoundError.")
            raise FileNotFoundError(f"The directory {outdir} does not exist. Please create it before running the script.")
        else:
            print(f"Directory exists: {outdir}")
        
        output_file = f'{outdir}/{c}.bed'
        print(f"Saving BED data to: {output_file}")
        np.savetxt(output_file, bed, fmt='%s', delimiter='\t')
        print(f"Saved BED data for chromosome {c}.")

print("Finished process_chromsizes.py script.")

