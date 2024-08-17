import gzip
import h5py
import numpy as np
import os
from config_and_print import filtered_list, resolutions, output_directory

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

prefix_file_path = filtered_list

def combine_tensors_for_chromosome(chromosome, prefix_file_path, resolution_label, output_directory):
    # Read prefixes from the file
    with open(prefix_file_path, 'r') as f:
        prefixes = [line.strip() for line in f]

    # Initialize a list to store the individual tensors
    tensors = []

    # Load each tensor in the order specified by the prefix list
    for prefix in prefixes:
        tensor_file_path = os.path.join(output_directory, f'hic_methy_{resolution_label}_tensor_singlecell/{chromosome}/{prefix}_{chromosome}.h5')
        
        if not os.path.exists(tensor_file_path):
            raise FileNotFoundError(f"Tensor file {tensor_file_path} not found.")
        
        with h5py.File(tensor_file_path, 'r') as hf:
            tensor = hf['Tensor'][:]
            tensors.append(tensor)

    # Stack the tensors along a new axis to create an n by n by 2 by m tensor
    combined_tensor = np.stack(tensors, axis=-1)

    # Define the output path for the combined tensor
    combined_tensor_path = os.path.join(output_directory, f'hic_methy_{resolution_label}_all_cells_tensors/{chromosome}_all_cells_tensor.h5')

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(combined_tensor_path), exist_ok=True)

    # Save the combined tensor to an HDF5 file
    with h5py.File(combined_tensor_path, 'w') as hf:
        hf.create_dataset('Combined_Tensor', data=combined_tensor)

    print(f"Combined tensor saved to {combined_tensor_path}")

# Example usage
for i in range(1, 23):  # Loop through chromosomes 1 to 22
    chromosome = f'chr{i}'
    combine_tensors_for_chromosome(chromosome, filtered_list, resolution_label, output_directory)













