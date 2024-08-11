#################################################################################
#create dark bins file if not already created
#################################################################################


############################################
# imports
############################################

import pyBigWig
import scipy.io as sio
import numpy as np
import math 
import matplotlib.pyplot as plt
import os
import pandas as pd
from heapq import nlargest
import copy
import matplotlib.gridspec as gridspec
import pandas as pd
import pickle
import seaborn as sns
import h5py
from scipy.stats import pearsonr
from config_and_print import methy_directory, filtered_list, chrom_file, resolutions, output_directory, mappability_threshold
#chromosomes = [f'chr{chrom}' for chrom in chromosomes]

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
    bigwig_file = "../../bin/softwarefiles/dark_regions_hg19.bigWig"
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
