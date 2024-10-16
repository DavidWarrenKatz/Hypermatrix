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

def get_matrix(file_name):
    # Construct the file path
    csr_mat = load_csr_matrix_from_hdf5(file_name)
    # Convert CSR matrix to a dense N-dimensional NumPy array
    dense_matrix = csr_mat.toarray()
    symmetric_matrix = dense_matrix + dense_matrix.T - np.diag(dense_matrix.diagonal())
    return symmetric_matrix

def emphasize_interactions(matrix, max_distance):
    """Emphasize interactions by adding offsets of contacts."""
    emphasized_matrix = csr_matrix(matrix.shape)
    for offset in range(1, max_distance + 1):
        emphasized_matrix += triu(matrix, offset) + tril(matrix, -offset)
    return emphasized_matrix

def csr_pearson_correlation(csr_mat):
    """Calculate Pearson correlation from a CSR matrix."""
    csc_mat = csr_mat.tocsc()
    
    mean = np.array(csc_mat.mean(axis=1)).flatten()
    std_dev = np.sqrt(csc_mat.power(2).mean(axis=1).A1 - mean**2)
    epsilon = 1e-10
    std_dev[std_dev == 0] = epsilon

    rows, cols = csr_mat.nonzero()
    standardized_data = (csr_mat.data - mean[rows]) / std_dev[rows]
    standardized_csr = csr_matrix((standardized_data, (rows, cols)), shape=csr_mat.shape)

    correlation_matrix = standardized_csr.dot(standardized_csr.T).toarray()
    diag = np.sqrt(np.diag(correlation_matrix))
    diag[diag == 0] = epsilon

    correlation_matrix /= diag[:, None]
    correlation_matrix /= diag[None, :]

    return csr_matrix(np.nan_to_num(correlation_matrix))


def process_matrices(input_dir, output_emphasized_dir, max_distance):
    """Process each Hi-C data file, compute matrices, and save the results."""
    os.makedirs(output_emphasized_dir, exist_ok=True)
    print(f"Input directory: {input_dir}")
    
    for file_path in glob.glob(os.path.join(input_dir, '*_pad1_std1_rp0.5_sqrtvc.hdf5')):
        print(f"Processing file: {file_path}")
        
        # Extract base name and modify it to desired format
        base_name = os.path.basename(file_path)
        
        # Split at the first occurrence of '_chr' to separate the prefix and the chromosome information
        prefix, chromosome_info = base_name.split('_chr', 1)
        
        # Further split the chromosome information to get the chromosome part only
        chromosome = 'chr' + chromosome_info.split('_')[0]  # Extracts chromosome like 'chr1'
        
        # Construct the desired output file name
        output_file_name = f"{prefix}_{chromosome}.h5"
        output_emphasized_path = os.path.join(output_emphasized_dir, output_file_name)
        
        # Skip computation if output file already exists
        if os.path.exists(output_emphasized_path):
            print(f"Skipping {output_emphasized_path}, emphasized Hi-C matrix already exists.")
            continue

        data = get_matrix(file_path)

        if not data.any():
            print(f"No data loaded from {file_path}")
            continue

        csr_mat = csr_matrix(data)
        if csr_mat is None:
            print(f"Failed to create matrix from data in {file_path}")
            continue
        
        emphasized_matrix = emphasize_interactions(csr_mat, max_distance)

        with h5py.File(output_emphasized_path, 'w') as output_file:
            grp = output_file.create_group('Matrix')
            grp.create_dataset('data', data=emphasized_matrix.data)
            grp.create_dataset('indices', data=emphasized_matrix.indices)
            grp.create_dataset('indptr', data=emphasized_matrix.indptr)
            grp.attrs['shape'] = emphasized_matrix.shape





max_genomic_distance = int(10_000_000 / resolution) + 1
base_input_dir = f'{output_directory}/hicluster_{resolution_label}_impute_dir/'
base_output_emphasized_dir = f'{output_directory}/hic_{resolution_label}_imputed_emphasized_dir/'

for i in range(1, 23):  
    chromosome = f'chr{i}'
    input_dir = base_input_dir + f'{chromosome}'
    output_emphasized_dir = base_output_emphasized_dir + f'{chromosome}'
    print(f'Processing {chromosome}')
    process_matrices(input_dir, output_emphasized_dir, max_genomic_distance)

