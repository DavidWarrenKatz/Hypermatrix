import h5py
import numpy as np
import os
from scipy.sparse import csr_matrix
from sklearn.decomposition import NMF
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

def normalize_matrix_columns(A):
    column_norms = np.linalg.norm(A, axis=0)
    zero_norms = column_norms == 0
    column_norms[zero_norms] = 1
    normalized_A = A / column_norms
    normalized_A[:, zero_norms] = 0
    return normalized_A

def load_csr_matrix_from_hdf5(file_path):
    with h5py.File(file_path, 'r') as file:
        data = file['Matrix/data'][:]
        indices = file['Matrix/indices'][:]
        indptr = file['Matrix/indptr'][:]
        shape = file['Matrix'].attrs['shape']
    return csr_matrix((data, indices, indptr), shape=shape)

def nmf(matrix, rank):
    model = NMF(n_components=rank, init='random', random_state=0)
    W = model.fit_transform(matrix)
    H = model.components_
    reconstructed_matrix = np.dot(W, H)
    return reconstructed_matrix

def normalize_matrix(A):
    norm = np.linalg.norm(A)
    if norm == 0:
        return A
    return A / norm

def process_tensor(methy_file_path, hic_path, output_tensor_path):
    with h5py.File(methy_file_path, 'r') as hf:
        dataset_name = list(hf.keys())[0]
        methy_matrix = hf[dataset_name][:]
    
    # Load the Hi-C matrix
    hic_matrix = load_csr_matrix_from_hdf5(hic_path)

    # Convert to dense matrix and drop the last row and column
    hic_matrix_dense = hic_matrix.toarray()
    hic_matrix_dense = hic_matrix_dense[:-1, :-1]

    # Update the methy_matrix size to match hic_matrix
    methy_matrix = methy_matrix[:hic_matrix_dense.shape[0], :hic_matrix_dense.shape[1]]

    # For all zero entries in hic_matrix, make corresponding entries in methy_matrix zero
    methy_matrix[hic_matrix_dense == 0] = 0

    # Normalize each matrix so that each column has a norm of 1
    hic_matrix = normalize_matrix_columns(hic_matrix_dense)
    methy_matrix = normalize_matrix_columns(methy_matrix)

    # This is my imputation step
    # Both matrices should be non-negative still, have not taken correlations
    # Take the NMF rank 5 of each matrix
    nmf_rank = 5
    hic_matrix_nmf = nmf(hic_matrix, nmf_rank)
    methy_matrix_nmf = nmf(methy_matrix, nmf_rank)

    #Take correlations and add by one to make non-negative
    hic_correlation = np.corrcoef(hic_matrix_nmf) + 1
    methy_correlation = np.corrcoef(methy_matrix_nmf) + 1

    # Additional normalization step
    # Both matrices must have the same norm before tensor decomposition
    hic_matrix_nmf = normalize_matrix(hic_correlation)
    methy_matrix_nmf = normalize_matrix(methy_correlation)

    # Combine both normalized matrices into a tensor
    tensor = np.stack((hic_matrix_nmf, methy_matrix_nmf), axis=-1)

    # Create the output directory if it does not exist
    os.makedirs(os.path.dirname(output_tensor_path), exist_ok=True)

    # Save the tensor to an HDF5 file
    with h5py.File(output_tensor_path, 'w') as hf:
        hf.create_dataset('Tensor', data=tensor)
    print(f"Tensor saved to {output_tensor_path}")

prefix_file_path = filtered_list
# Read prefixes from the file
with open(prefix_file_path, 'r') as f:
    prefixes = [line.strip() for line in f]    

for i in range(1, 23):  
    chromosome = f'chr{i}'
    for prefix in prefixes:
        methy_file_path = f'{output_directory}/methy_{resolution_label}_outerproduct_dir/{chromosome}/{prefix}_outer_product.h5'
        hic_path = f'{output_directory}/hic_{resolution_label}_emphasized_dir/{chromosome}/{prefix}_{chromosome}.h5'
        output_tensor_basepath = f'{output_directory}/hic_methy_{resolution_label}_tensor_singlecell/{chromosome}/'
        output_tensor_path = output_tensor_basepath + f'{prefix}_{chromosome}.h5'
        
        if os.path.exists(output_tensor_path):
            print(f"File {output_tensor_path} already exists. Skipping.")
        else:
            # Ensure the output directory exists
            os.makedirs(output_tensor_basepath, exist_ok=True)
            process_tensor(methy_file_path, hic_path, output_tensor_path)

