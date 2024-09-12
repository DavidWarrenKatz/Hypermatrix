#################################################################################
#create dark bins file if not already created
#################################################################################


############################################
# imports
############################################

import math
import numpy as np
import pyBigWig

import sys
import os

# Add the directory where config.py is located to the Python path
config_dir = '../../../'
config_dir = os.path.abspath(config_dir)  # Get absolute path

config_file = os.path.join(config_dir, 'config.py')  # Full path to config.py

# Check if the directory and config.py exist
if os.path.isdir(config_dir) and os.path.isfile(config_file):
    sys.path.append(config_dir)
    print(f"Config directory added to sys.path: {config_dir}")
    print(f"Found config.py at: {config_file}")
else:
    raise FileNotFoundError(f"config.py not found in directory: {config_dir}")

from config import methy_directory, filtered_list, chrom_file, resolutions, output_directory, reference_genome, mappability_threshold

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

# Check if the bins to remove file has already been created
bins_file_path = f'{output_directory}/bins_to_remove_res{resolution_label}.npz'
if os.path.exists(bins_file_path):
    print(f"{bins_file_path} already exists. Skipping computation.")	
else:
    bigwig_file = "{output_directory}/../../src/softwarefiles/dark_regions_hg19.bigWig"
    # Open the BigWig file
    bw = pyBigWig.open(bigwig_file)

    # Define the chromosomes you want to analyze
    chromosomes = ['chr' + str(i) for i in range(1, 23)] 

    # Define the threshold for removing bins based on average mappability
    threshold = mappability_threshold

    # Create a dictionary to store the bin indices to remove for each chromosome
    bins_to_remove = {}

    # Loop through each chromosome
    for chrom in chromosomes:
        chrom_size = bw.chroms(chrom)

        if chrom_size is None:
            print(f"Chromosome {chrom} not found in the BigWig file.")
            continue

        # Calculate the number of bins based on the specified resolution
        num_bins = math.ceil(chrom_size / resolution) #last bin may not be of size resolution

        # Create lists to store bin indices to remove
        remove_indices = []

        # Calculate average mappability for each bin
        for i in range(num_bins):
            # Determine the start and end positions of the bin
            start = i * resolution
            end = min((i + 1) * resolution, chrom_size)  # to account for last bin which may be incomplete

            # Extract the mappability values for the bin
            values = np.nan_to_num(bw.values(chrom, start, end))

            # Calculate the average mappability score for the bin
            avg_mappability = np.mean(values)

            # Check if the average mappability is below the threshold
            if avg_mappability < threshold:
                remove_indices.append(i)

        # Store the bin indices to remove for this chromosome
        bins_to_remove[chrom] = remove_indices

    # Close the BigWig file
    bw.close()

    # Convert the lists in bins_to_remove to numpy arrays
    for chrom in bins_to_remove:
        bins_to_remove[chrom] = np.array(bins_to_remove[chrom])

    # Save the dictionary as an .npz file
    np.savez(bins_file_path, **bins_to_remove)
    print(f"Bins to remove file created and saved to {bins_file_path}")
