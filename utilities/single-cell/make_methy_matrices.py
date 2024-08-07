import gzip
import h5py
import numpy as np
import os
from config_and_print import methy_directory, filtered_list, chrom_file, resolutions, output_directory

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]
resolution_value, resolution_label = resolution_str.split(':')

# Convert resolution value to integer
resolution = int(resolution_value)

methy_output_dir = os.path.join(output_directory, f"methy_{resolution_label}_outerproduct_dir")
methy_matrix_path = os.path.join(output_directory, f'b37.autosome.{resolution_label}_interval.add_value.methy.bed.gz')
prefix_file_path = filtered_list

# Function to normalize each column of a matrix so that each column has a norm of one.
def normalize_matrix_columns(A):
    column_norms = np.linalg.norm(A, axis=0)
    if np.any(column_norms == 0):
        raise ValueError("One or more columns have zero norm. Cannot normalize those columns.")
    normalized_A = A / column_norms
    return normalized_A

# Function to calculate the matrix from a compressed file.
def calculate_matrix(file_path):
    data = []
    with gzip.open(file_path, 'rt') as file:
        for line in file:
            fields = line.strip().split('\t')
            numerator_indices = range(6, len(fields), 2)
            denominator_indices = range(7, len(fields), 2)
            row_data = []
            for num_idx, denom_idx in zip(numerator_indices, denominator_indices):
                numerator = float(fields[num_idx])
                denominator = float(fields[denom_idx])
                row_data.append(0 if denominator == 0 else numerator / denominator)
            data.append(row_data)
    data_matrix = np.array(data)
    return data_matrix

# Read chromosome sizes from the chrom_file
chromosome_lengths = {}
with open(chrom_file, 'r') as f:
    for line in f:
        chrom, length = line.strip().split()
        chromosome_lengths[chrom] = int(length)

# Normalize the methylation matrix columns
methylation_matrix = normalize_matrix_columns(calculate_matrix(methy_matrix_path))
print(f'methylation matrix shape at resolution {resolution}')
print(methylation_matrix.shape)

# Calculate bins per chromosome and slice data accordingly
chromosome_bins = {key: length // resolution for key, length in chromosome_lengths.items()}
chromosome_data = {}
start_bin = 0

for chromosome, bins in chromosome_bins.items():
    end_bin = start_bin + bins
    chromosome_data[chromosome] = methylation_matrix[start_bin:end_bin, :]
    start_bin = end_bin

# Ensure the output directory exists
os.makedirs(methy_output_dir, exist_ok=True)

# Read prefixes from the file
with open(prefix_file_path, 'r') as f:
    prefixes = [line.strip() for line in f]

# File to save all tensors
hdf5_filename = os.path.join(methy_output_dir, f"all_chromosomes_methylation_tensor_{resolution_label}.h5")

if not os.path.exists(hdf5_filename):
    with h5py.File(hdf5_filename, 'w') as all_chromosomes_file:
        # Loop through each chromosome and process its data
        for chromosome, data in chromosome_data.items():
            num_samples = data.shape[1]  # Number of samples
            outer_product_matrices = []

            # Create directory for the chromosome
            chromosome_dir = os.path.join(methy_output_dir, f"chr{chromosome}")
            os.makedirs(chromosome_dir, exist_ok=True)

            # Compute the outer product for each sample
            for sample_index in range(num_samples):
                vector = data[:, sample_index]
                outer_product_matrix = np.outer(vector, vector)
                outer_product_matrices.append(outer_product_matrix)

                # Save each outer product matrix in its own .h5 file if it doesn't exist
                prefix = prefixes[sample_index]
                individual_matrix_filename = os.path.join(chromosome_dir, f"{prefix}_outer_product.h5")
                if not os.path.exists(individual_matrix_filename):
                    with h5py.File(individual_matrix_filename, 'w') as individual_matrix_file:
                        individual_matrix_file.create_dataset(f"chr_{chromosome}", data=outer_product_matrix, compression="gzip")

            # Stack all the matrices into a tensor
            tensor = np.stack(outer_product_matrices, axis=0)

            # Transpose the tensor to put the sample index as the third dimension
            tensor_transposed = tensor.transpose((1, 2, 0))

            # Save the transposed tensor to the collective HDF5 file
            all_chromosomes_file.create_dataset(f"chr{chromosome}", data=tensor_transposed, compression="gzip")

            # Also save to an individual HDF5 file for this chromosome if it doesn't exist
            individual_hdf5_filename = os.path.join(chromosome_dir, f"chr{chromosome}_tensor_methylation.h5")
            if not os.path.exists(individual_hdf5_filename):
                with h5py.File(individual_hdf5_filename, 'w') as individual_chromosome_file:
                    individual_chromosome_file.create_dataset(f"chr_{chromosome}", data=tensor_transposed, compression="gzip")

            print(f"Saved tensor for chromosome {chromosome} in both collective and individual HDF5 files.")
else:
    print(f"File {hdf5_filename} already exists. Skipping computation.")

