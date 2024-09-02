import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.sparse import csr_matrix, triu, tril, csc_matrix
import os
import glob
from config_and_print import methy_directory, filtered_list, chrom_file, resolutions, output_directory

# Ensure resolutions is treated as a tuple or list of strings
if isinstance(resolutions, str):
    resolutions = (resolutions,)

# Print resolutions for debugging
print(f"Resolutions from config: {resolutions}")

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

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]
resolution, resolution_label = parse_resolution(resolution_str)

print(f"Extracted resolution: {resolution}, label: {resolution_label}")

def load_csr_matrix_from_hdf5(file_path):
    """Load a CSR matrix from an HDF5 file."""
    with h5py.File(file_path, 'r') as file:
        data = file['Matrix/data'][:]
        indices = file['Matrix/indices'][:]
        indptr = file['Matrix/indptr'][:]
        shape = file['Matrix'].attrs['shape']
    return csr_matrix((data, indices, indptr), shape=shape)

def third_order_cumulant_matrix(data):
    """Compute third-order cumulant tensor."""
    if not isinstance(data, np.ndarray):
        data = data.toarray()

    symmetric_matrix = data + data.T - np.diag(data.diagonal())
    n_columns = symmetric_matrix.shape[1]
    means = np.mean(symmetric_matrix, axis=0)
    
    cumulants = np.zeros((n_columns, n_columns, n_columns))
    
    for i in range(n_columns):
        for j in range(i, n_columns):
            for k in range(j, n_columns):
                x, y, z = symmetric_matrix[:, i], symmetric_matrix[:, j], symmetric_matrix[:, k]
                cumulant_ijk = np.mean((x - means[i]) * (y - means[j]) * (z - means[k])) - \
                               means[i] * np.mean((y - means[j]) * (z - means[k])) - \
                               means[j] * np.mean((x - means[i]) * (z - means[k])) - \
                               means[k] * np.mean((x - means[i]) * (y - means[j])) + \
                               2 * means[i] * means[j] * means[k]

                cumulants[i, j, k] = cumulants[i, k, j] = cumulants[j, i, k] = \
                                     cumulants[j, k, i] = cumulants[k, i, j] = cumulants[k, j, i] = cumulant_ijk
    
    return cumulants

def process_hic_files_for_cumulant(input_dir, output_cumulant_dir):
    """Process each Hi-C file and calculate the third-order cumulant tensor."""
    os.makedirs(output_cumulant_dir, exist_ok=True)
    
    for file_path in glob.glob(os.path.join(input_dir, '*.h5')):
        file_name = os.path.splitext(os.path.basename(file_path))[0] + '_cumulant.h5'
        output_path = os.path.join(output_cumulant_dir, file_name)
        
        # Check if the cumulant tensor already exists
        if os.path.exists(output_path):
            print(f"Cumulant tensor already exists for {file_name}. Skipping computation.")
            continue
        
        print(f"Processing Hi-C matrix from file: {file_path}")
        
        csr_mat = load_csr_matrix_from_hdf5(file_path)
        
        if csr_mat is None:
            print(f"Failed to load CSR matrix from {file_path}")
            continue
        
        cumulant_tensor = third_order_cumulant_matrix(csr_mat)
        print(f"Cumulant tensor for {file_path} calculated.")
        
        with h5py.File(output_path, 'w') as h5file:
            h5file.create_dataset('cumulant_tensor', data=cumulant_tensor)
        print(f"Cumulant tensor saved to {output_path}")

base_output_emphasized_dir = f'{output_directory}/hic_{resolution_label}_emphasized_dir/'
base_output_cumulant_dir = f'{output_directory}/hic_{resolution_label}_cumulant_dir/'

for i in range(1, 23):  
    chromosome = f'chr{i}'
    input_dir = os.path.join(base_output_emphasized_dir, chromosome)
    output_cumulant_dir = os.path.join(base_output_cumulant_dir, chromosome)
    
    print(f'Processing cumulant tensors for {chromosome}')
    process_hic_files_for_cumulant(input_dir, output_cumulant_dir)

