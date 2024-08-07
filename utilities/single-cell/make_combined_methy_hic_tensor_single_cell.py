import h5py
import gzip
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from config_and_print import methy_directory, filtered_list, chrom_file, resolutions, output_directory

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

refix_file_path = filtered_list

# File to save all tensors
hdf5_filename = os.path.join(methy_output_dir, f"all_chromosomes_methylation_tensor_{resolution_label}.h5")

# Check if the file already exists and skip computations if it does
if os.path.exists(hdf5_filename):
    print(f"File {hdf5_filename} already exists. Skipping all computations.")
else:
    # Ensure the output directory exists
    os.makedirs(methy_output_dir, exist_ok=True)

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
    print(f'methylation matrix shape at resolution {resolution_label}')
    print(methylation_matrix.shape)

    # Calculate bins per chromosome and slice data accordingly
    chromosome_bins = {key: length // resolution for key, length in chromosome_lengths.items()}
    chromosome_data = {}
    start_bin = 0

    for chromosome, bins in chromosome_bins.items():
        end_bin = start_bin + bins
        chromosome_data[chromosome] = methylation_matrix[start_bin:end_bin, :]
        start_bin = end_bin

    # Read prefixes from the file
    with open(prefix_file_path, 'r') as f:
        prefixes = [line.strip() for line in f]

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

# Function to normalize each column of a matrix so that each column has a norm of one.
def normalize_matrix_columns(A):
    column_norms = np.linalg.norm(A, axis=0)
    # Avoid division by zero by setting zero norms to 1 temporarily
    zero_norms = column_norms == 0
    column_norms[zero_norms] = 1
    normalized_A = A / column_norms
    # Restore zero columns to zeros
    normalized_A[:, zero_norms] = 0
    return normalized_A

def load_csr_matrix_from_hdf5(file_path):
    with h5py.File(file_path, 'r') as file:
        data = file['Matrix/data'][:]
        indices = file['Matrix/indices'][:]
        indptr = file['Matrix/indptr'][:]
        shape = file['Matrix'].attrs['shape']  # Assuming shape is stored as an attribute
    return csr_matrix((data, indices, indptr), shape=shape)

# Read chromosome sizes from the chrom_file
chromosome_lengths = {}
with open(chrom_file, 'r') as f:
    for line in f:
        chrom, length = line.strip().split()
        chromosome_lengths[chrom] = int(length)

chromosome_bins = {key: length // resolution for key, length in chromosome_lengths.items()}

# Read prefixes from the file
with open(prefix_file_path, 'r') as f:
    prefixes = [line.strip() for line in f]

for chromosome, bins in chromosome_bins.items():
    for prefix in prefixes:
        methy_file_path = f'{output_directory}/methy_{resolution_label}_outerproduct_dir/chr{chromosome}/{prefix}_outer_product.h5'
        hic_path = f'{output_directory}/hic_{resolution_label}_correlation_dir/chr{chromosome}/{prefix}_chr{chromosome}.h5'

        if not os.path.exists(methy_file_path) or not os.path.exists(hic_path):
            print(f"Files {methy_file_path} or {hic_path} do not exist. Skipping.")
            continue

        with h5py.File(methy_file_path, 'r') as hf:
            dataset_name = list(hf.keys())[0]  
            methy_matrix = hf[dataset_name][:]
        
        # Check for any negative values in the matrix
        if np.any(methy_matrix < 0):
            print(f"Warning: The matrix contains negative values.")

        # Load the Hi-C matrix
        hic_matrix = load_csr_matrix_from_hdf5(hic_path)

        # Convert to dense matrix and drop the last row and column
        hic_matrix_dense = hic_matrix.toarray()
        hic_matrix_dense = hic_matrix_dense[:-1, :-1]

        # Update the methy_matrix size to match hic_matrix
        methy_matrix = methy_matrix[:hic_matrix_dense.shape[0], :hic_matrix_dense.shape[1]]

        # For all zero entries in hic_matrix, make corresponding entries in methy_matrix zero
        methy_matrix[hic_matrix_dense == 0] = 0

        # Normalize each matrix so that each matrix has a norm of 1
        hic_matrix_dense = normalize_matrix_columns(hic_matrix_dense)
        methy_matrix = normalize_matrix_columns(methy_matrix)

        # Combine both normalized matrices into a tensor
        tensor = np.stack((hic_matrix_dense, methy_matrix), axis=-1)

        # Create the directory if it does not exist
        output_tensor_path = f'{output_directory}/combined_{resolution_label}_methy_hic_single_cell_tensor/chr{chromosome}/{prefix}.h5'
        os.makedirs(os.path.dirname(output_tensor_path), exist_ok=True)

        # Save the tensor to an HDF5 file
        with h5py.File(output_tensor_path, 'w') as hf:
            hf.create_dataset('Tensor', data=tensor)

        print(f"Tensor saved to {output_tensor_path}")

        # Plotting the slices of the tensor
        plt.figure(figsize=(10, 5))

        plt.subplot(1, 2, 1)
        plt.imshow(tensor[:, :, 0], aspect='auto', cmap='viridis')
        plt.colorbar()
        plt.title(f'Hi-C Matrix Slice - {chromosome} - {prefix}')

        plt.subplot(1, 2, 2)
        plt.imshow(tensor[:, :, 1], aspect='auto', cmap='viridis')
        plt.colorbar()
        plt.title(f'Methylation Matrix Slice - {chromosome} - {prefix}')

        plt.show()

        # Calculating and printing statistics
        def matrix_statistics(matrix, name):
            shape = matrix.shape
            min_val = np.min(matrix)
            max_val = np.max(matrix)
            num_nonzeros = np.count_nonzero(matrix)
            norm = np.linalg.norm(matrix)
            column_norms = np.linalg.norm(matrix, axis=0)

            print(f"Statistics for {name}:")
            print(f"Shape: {shape}")
            print(f"Min: {min_val}")
            print(f"Max: {max_val}")
            print(f"Number of non-zeros: {num_nonzeros}")
            print(f"Norm: {norm}")
            print(f"Norm of each column: {column_norms}")
            print()

        matrix_statistics(hic_matrix_dense, f"Hi-C Matrix - {chromosome} - {prefix}")
        matrix_statistics(methy_matrix, f"Methylation Matrix - {chromosome} - {prefix}")

