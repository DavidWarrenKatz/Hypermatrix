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
from config_and_print import methy_directory, filtered_list, chrom_file, resolutions, output_directory, mappability_threshold, normalization
#chromosomes = [f'chr{chrom}' for chrom in chromosomes]

# Ensure resolutions is treated as a tuple or list of strings
if isinstance(resolutions, str):
    resolutions = (resolutions,)

# Print resolutions for debugging
print(f"Resolutions from config: {resolutions}")

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]

# Debug print to check the value of resolution_str
print(f"Extracted resolution string: {resolution_str}")

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

########################################################################
# create the cell type dictionary
# [TO DO] This needs to be replaced with SNPS code 
########################################################################
# Define the path file with prefixes and colors in the following form
#1       sc1.ACTTGA      red
#2       sc1.GCCAAT      red
#3       sc1.TAGCTT      red
#4       sc10.TAGCTT     blue
#
filename = '../../bin/name.order.HCG_methy.with_color.txt'

# Initialize an empty dictionary to store cell ID and color
cell_color_dict = {}

# Open and read the file
with open(filename, 'r') as file:
    for line in file:
        # Split the line into parts
        parts = line.strip().split()
        # Extract cell ID and color
        cell_id = parts[1]
        color = parts[2]
        # Store in dictionary
        cell_color_dict[cell_id] = color

# Define the path to the tensor sample order file
#This file contains the prefixes in the form
#sc11.ACTTGA
#sc11.CGATGT
#sc11.GCCAAT
#
tensor_order_filename = f'{output_directory}/filtered_bam_list.txt'

# Initialize a list to store the 1s and 0s
color_vector = []

# Open and read the tensor sample order file
with open(tensor_order_filename, 'r') as file:
    for line in file:
        sample_id = line.strip()  # Remove any trailing newlines or spaces
        if sample_id in cell_color_dict and cell_color_dict[sample_id] == 'red':
            color_vector.append(1)
        else:
            color_vector.append(0)

# Output the color vector to check
print(len(color_vector))

# Create a mapping dictionary
color_mapping = {
    'red': 'imr90',
    'blue': 'gm12878'
}

# Update the dictionary using the mapping
updated_cell_color_dict = {key: color_mapping[value] for key, value in cell_color_dict.items()}

##############################################################################
#import bins to remove created previously
#############################################################################

#Import bins to remove dictionary
bins_file_path = f'{output_directory}/bins_to_remove_res{resolution_label}.npz'
if os.path.exists(bins_file_path):
    print(f"{bins_file_path} exist. Importing.")
    #load the dark regions data and the A/B compartment data
    loaded_data = np.load(bins_file_path)
    # Convert the loaded data back into a dictionary with the same structure
    bins_to_remove = {chrom: loaded_data[chrom] for chrom in loaded_data}

###########################################################################    
#create a dictionary of the A/B compartment calls for the bulk data
###########################################################################
bulk_data = {}
path_to_eigenvectors = '../../projects/single_cell_files/eigenvector/'

for i in range(1, 23):
    file = path_to_eigenvectors + f'res{resolution}_ch{i}_oe_GM12878_{normalization}_eigenvector.txt'
    key = os.path.splitext(os.path.basename(file))[0]  
    bulk_data[key] = pd.read_csv(file, header=None, names=['eigenvalue'])
    file = path_to_eigenvectors + f'res{resolution}_ch{i}_oe_IMR90_{normalization}_eigenvector.txt'
    key = os.path.splitext(os.path.basename(file))[0]  
    bulk_data[key] = pd.read_csv(file, header=None, names=['eigenvalue'])
    
###########################################################################
# download H3K9ac file if it does not exist
###########################################################################

def calculate_bin_averages(data, elements_per_bin):
    num_bins = math.ceil(len(data)/elements_per_bin)
    bin_averages = np.zeros(num_bins)
    for i in range(num_bins):
        start_index = i * elements_per_bin
        end_index = min((i + 1) * elements_per_bin, len(data))
        bin_data = data[start_index:end_index]
        if len(bin_data) > 0:
            bin_averages[i] = np.mean(bin_data)
        else:
            bin_averages[i] = 0
    return np.nan_to_num(bin_averages, nan=0.0)

# Read chromosome sizes from hg19.autosome.chrom.sizes
lengths = {}
with open(chrom_file, 'r') as file:
    for line in file:
        chrom, size = line.strip().split()
        lengths[chrom] = int(size)

# Check if the output file already exists
output_file = f'../../bin/softwarefiles/h3k9ac_res{resolution}_GM12878.pkl'

if not os.path.exists(output_file):
    # Initialization
    chromosomes = list(lengths.keys())
    h3k9ac = {name: [] for name in chromosomes}

    bigwig_H3K9ac_path = '../../bin/softwarefiles/ENCFF128UVW_hg19_H3K9ac_GM12878.bigWig'
    if not os.path.exists(bigwig_H3K9ac_path):
        os.system(f'wget https://www.encodeproject.org/files/ENCFF128UVW/@@download/ENCFF128UVW.bigWig -O {bigwig_H3K9ac_path}')

    bigwig_H3K9ac = pyBigWig.open(bigwig_H3K9ac_path)
    
    # Process each chromosome
    for chromosome_name in chromosomes:
        start_position = 1
        end_position = lengths[chromosome_name]
        values_file1 = np.array(bigwig_H3K9ac.values(chromosome_name, start_position, end_position))
        h3k9ac[chromosome_name] = calculate_bin_averages(values_file1, resolution)

    # Save the results in a pickle file
    with open(output_file, 'wb') as file:
        pickle.dump(h3k9ac, file)
else:
    print(f"The file {output_file} already exists.")
    with open(output_file, 'rb') as file:
        h3k9ac = pickle.load(file)

################################################################################    
#make sure each GM12878 eigenvector has positive value for active A compartment
################################################################################
for i in range(1, 23):
    key_gm12878 = f'res{resolution}_ch{i}_oe_GM12878_{normalization}_eigenvector'
    chromosome_key = f'chr{i}'
    h3k9ac_df = pd.DataFrame(h3k9ac[chromosome_key], columns=['H3K9ac_signal'])
    
    df_gm12878_positive = bulk_data[key_gm12878]['eigenvalue']  
    df_gm12878_negative = -bulk_data[key_gm12878]['eigenvalue']

    # Compute correlations by first ensuring eigenvector data is in DataFrame format
    corr_positive = df_gm12878_positive.corr(h3k9ac_df['H3K9ac_signal'])
    corr_negative = df_gm12878_negative.corr(h3k9ac_df['H3K9ac_signal'])
    if corr_negative > corr_positive:
        bulk_data[key_gm12878]['eigenvalue'] = -bulk_data[key_gm12878]['eigenvalue']
        print(f"Switched for chromosome {i}")    

###############################################################################    
#make sure each IMR90 eigenvector has consistent orientation with GM12878
################################################################################
for i in range(1, 23):
    # Construct keys for GM12878 and IMR90
    key_gm12878 = f'res{resolution}_ch{i}_oe_GM12878_{normalization}_eigenvector'
    key_imr90 = f'res{resolution}_ch{i}_oe_GM12878_{normalization}_eigenvector'
    
    # Retrieve DataFrames for GM12878 and IMR90
    df_gm12878 = bulk_data[key_gm12878]
    df_imr90 = bulk_data[key_imr90]
    
    # Ensure data is in expected format
    if not df_gm12878.empty and not df_imr90.empty:
        # Concatenate DataFrames side-by-side
        combined_df = pd.concat([df_gm12878.reset_index(drop=True), df_imr90.reset_index(drop=True)], axis=1, keys=['gm12878', 'imr90'])
        
        # Drop rows with NaN values in either column
        combined_df.dropna(inplace=True)
        
        # Extract Series after dropping NaNs
        gm12878_series = combined_df['gm12878']['eigenvalue']
        imr90_series = combined_df['imr90']['eigenvalue']

        # Calculate correlations
        corr_positive = gm12878_series.corr(imr90_series)
        corr_negative = gm12878_series.corr(-imr90_series)

        # If negating IMR90 improves correlation, update the original data
        if corr_negative > corr_positive:
            bulk_data[key_imr90]['eigenvalue'] = -bulk_data[key_imr90]['eigenvalue']
            print(f"switched for chromosome {i}")

#The A/B compartments of the proper orientation before dropping dark regions
original_bulk_data = copy.deepcopy(bulk_data)

#################################################################################
#remove dark regions
#dark reigons are obviously correlated
#I want to remove dark regions to get meeaningfully correlated regions
###############################################################################

for i in range(1, 23):
    # Construct the keys for GM12878 and IMR90
    key_gm12878 = f'res{resolution}_ch{i}_oe_GM12878_{normalization}_eigenvector'
    key_imr90 = f'res{resolution}_ch{i}_oe_IMR90_{normalization}_eigenvector'
    chrom = f'chr{i}'

    # Check if the chromosome exists in the bins_to_remove and in the data
    if chrom in bins_to_remove and key_gm12878 in bulk_data and key_imr90 in bulk_data:
        # Retrieve the indices to remove for this chromosome
        indices_to_remove = bins_to_remove[chrom]

        # Initialize lists to hold valid indices for both DataFrames
        valid_indices = []
        
        for idx in indices_to_remove:
            if idx < len(bulk_data[key_gm12878]) and idx < len(bulk_data[key_imr90]):
                valid_indices.append(idx)
        
        # Now drop the valid indices from both DataFrames
        if valid_indices:
            bulk_data[key_gm12878] = bulk_data[key_gm12878].drop(valid_indices).reset_index(drop=True)
            bulk_data[key_imr90] = bulk_data[key_imr90].drop(valid_indices).reset_index(drop=True)


# This is where AB_calls are saved

def extract_prefixes(file_path):
    """Extract prefixes from the filtered_bam_list.txt file."""
    prefixes = []
    with open(file_path, 'r') as file:
        for line in file:
            # Assuming the prefix format is like "sc11.ACTTGA" or similar
            parts = line.strip().split('.')
            if len(parts) >= 2:
                prefix = parts[0] + '.' + parts[1]  # Concatenate the first two parts to get the full prefix
                prefixes.append(prefix)
    return prefixes

def load_h5_file(file_path, dataset_name):
    """Load a dataset from an HDF5 file."""
    with h5py.File(file_path, 'r') as f:
        data = f[dataset_name][:]
    return data

def get_best_correlated_vector(V, eigenvector):
    """Find the row in V that best correlates with the given eigenvector."""
    eigenvector = eigenvector.values.flatten()
    best_corr = -np.inf  # Initialize with a very low correlation
    best_index = -1
    best_vector = None

    # Check each row in V for correlation with the eigenvector
    for i in range(V.shape[0]):  # Iterate over all rows in V
        # Synchronize non-NaN data
        valid_indices = ~np.isnan(V[i, :]) & ~np.isnan(eigenvector)
        if np.any(valid_indices):
            corr, _ = pearsonr(V[i, valid_indices], eigenvector[valid_indices])
            if corr > best_corr:  # Check if this is the best correlation so far
                best_corr = corr
                best_index = i
                best_vector = V[i, :]

    return best_vector, best_index, best_corr

# Directory setup
output_directory = '../../projects/single_cell_files/'
filtered_bam_list = os.path.join(output_directory, 'filtered_bam_list.txt')
base_tensor_dir = os.path.join(output_directory, 'tensor_1Mb_AB_factors')

# Extract prefixes from the file
prefixes = extract_prefixes(filtered_bam_list)
print(f"Extracted prefixes: {prefixes}")

# Dictionary to hold the dataframes for each chromosome
chromosome_results = {}

def normalize_vectors(V):
    """Normalize the columns of V to have a norm of 1."""
    norms = np.linalg.norm(V, axis=0)
    norms[norms == 0] = 1  # Prevent division by zero
    V_normalized = V / norms
    return .02 * V_normalized - .01

# Function to replace specific bins with NaNs, ensuring the bin is valid
def replace_bins_with_nans(compartment_values, bins_to_remove):
    compartment_values = np.array(compartment_values)  # Convert to numpy array if not already
    valid_bins_to_remove = [bin_idx for bin_idx in bins_to_remove if bin_idx < len(compartment_values)]
    compartment_values[valid_bins_to_remove] = np.nan  # Replace values in valid bins to remove with NaNs
    return compartment_values

# Loop through each chromosome directory
for i in range(1, 23):
    chromosome = f'chr{i}'
    results_df = pd.DataFrame(columns=['Sample', 'A/B Compartment', 'Correlation With Bulk', 'Cell Type'])

    for prefix in prefixes:
        input_file = os.path.join(base_tensor_dir, chromosome, f'{prefix}_compartments.h5')
        output_dir = os.path.join(output_directory, 'tensor_1Mb_AB_calls', chromosome)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_file = os.path.join(output_dir, f'{prefix}_tensor_AB_compartment_call.h5')
        text_output_file = os.path.join(output_dir, f'{prefix}_tensor_AB_compartment_call.txt')

        if os.path.exists(output_file):
            print(f"{output_file} already exists. Loading existing results.")
            with h5py.File(output_file, 'r') as output_h5:
                best_vector = output_h5['AB_Compartment'][:]
                best_corr = output_h5['correlation']
                cell_type = output_h5['cell_type']
 
        else:
            sample_id = prefix  
            cell_type = updated_cell_color_dict.get(sample_id, 'GM12878')  # Default set to GM12878
            key = f'res{resolution}_ch{i}_oe_{cell_type.upper()}_{normalization}_eigenvector'  

            print(f"Processing sample_id: {sample_id}, cell_type: {cell_type}, key: {key}, resolution {resolution}")

            if key in original_bulk_data:
                tensor_factors = load_h5_file(input_file, '/compartment_factors')
                tensor_factors = normalize_vectors(tensor_factors)
                print(f"Loaded tensor factors from {input_file}, shape: {tensor_factors.shape}")
                bulk_eigenvector = original_bulk_data[key]['eigenvalue'][:-1]
                print(f"Bulk eigenvector shape: {bulk_eigenvector.shape}")

                best_vector, best_index, best_corr = get_best_correlated_vector(tensor_factors, bulk_eigenvector)
                print(f"Best correlation: {best_corr} at index {best_index}")

                # Retrieve the bins to remove for this chromosome
                bins_to_remove_for_chrom = bins_to_remove[chromosome]
        
                # Apply the remove dark bins function to each cell's A/B Compartment data
                best_vector = replace_bins_with_nans(best_vector, bins_to_remove_for_chrom)
                
                # Save the results to the output file (HDF5)
                with h5py.File(output_file, 'w') as output_h5:
                    output_h5.create_dataset('AB_Compartment', data=best_vector)
                    output_h5.create_dataset('correlation', data=best_corr)
                    output_h5.create_dataset('cell_type', data=cell_type)

                # Save the A/B Compartment values as a text file
                np.savetxt(text_output_file, best_vector, fmt='%f')
                print(f"Saved A/B Compartment values to {text_output_file}")

        new_row = pd.DataFrame({
            'Sample': [sample_id],
            'A/B Compartment': [best_vector],
            'Correlation With Bulk': [best_corr],
            'Cell Type': [cell_type]
        })
        results_df = pd.concat([results_df, new_row], ignore_index=True)

    print(f"Finished processing {chromosome}, {results_df.shape[0]} rows added.")
    chromosome_results[chromosome] = results_df
                                                    
# Save original_bulk_data to a file
bulk_data_output_file = '../../projects/single_cell_files/original_bulk_data.pkl'
with open(bulk_data_output_file, 'wb') as f:
    pickle.dump(original_bulk_data, f)
print(f"original_bulk_data saved to {bulk_data_output_file}")

# Save chromosome_results to a file
chromosome_results_output_file = '../../projects/single_cell_files/chromosome_results.pkl'
with open(chromosome_results_output_file, 'wb') as f:
    pickle.dump(chromosome_results, f)
print(f"chromosome_results saved to {chromosome_results_output_file}")                                                               
