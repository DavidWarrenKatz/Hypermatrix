import os
import glob
import h5py
import numpy as np
from scipy.sparse import csr_matrix
from config_and_print import filtered_list, chrom_file, resolutions, output_directory

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

resolution_str = resolutions[0]
resolution, resolution_label = parse_resolution(resolution_str)

print(f"Extracted resolution: {resolution}, label: {resolution_label}")

def load_matrix_from_hdf5(file_path):
    """Load a matrix from an HDF5 file."""
    try:
        with h5py.File(file_path, 'r') as file:
            matrix = file[list(file.keys())[0]][:]
        return matrix
    except Exception as e:
        print(f"Error loading matrix from {file_path}: {e}")
        return None

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

def process_methylation_files_for_cumulant(input_dir, output_cumulant_dir):
    """Process each methylation matrix file and calculate the third-order cumulant tensor."""
    os.makedirs(output_cumulant_dir, exist_ok=True)
    
    h5_files = glob.glob(os.path.join(input_dir, '*.h5'))
    if not h5_files:
        print(f"No HDF5 files found in directory: {input_dir}")
        return
    
    for file_path in h5_files:
        print(f"Processing methylation matrix from file: {file_path}")
        
        matrix = load_matrix_from_hdf5(file_path)
        
        if matrix is None:
            print(f"Failed to load matrix from {file_path}")
            continue
        
        print(f"Loaded matrix shape: {matrix.shape}")
        
        cumulant_tensor = third_order_cumulant_matrix(matrix)
        print(f"Cumulant tensor for {file_path} calculated.")
        
        file_name = os.path.splitext(os.path.basename(file_path))[0] + '_cumulant.h5'
        output_path = os.path.join(output_cumulant_dir, file_name)
        
        with h5py.File(output_path, 'w') as h5file:
            h5file.create_dataset('cumulant_tensor', data=cumulant_tensor)
        print(f"Cumulant tensor saved to {output_path}")

base_output_methy_cumulant_dir = f'{output_directory}/methy_{resolution_label}_cumulant_dir/'

for i in range(1, 23):  
    chromosome = f'chr{i}'
    input_dir = os.path.join(output_directory, f'methy_{resolution_label}_outerproduct_dir/{chromosome}')
    output_cumulant_dir = os.path.join(base_output_methy_cumulant_dir, chromosome)
    
    print(f'Processing cumulant tensors for methylation data on {chromosome}')
    process_methylation_files_for_cumulant(input_dir, output_cumulant_dir)

