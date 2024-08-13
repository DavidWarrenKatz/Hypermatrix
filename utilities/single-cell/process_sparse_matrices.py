import os
import numpy as np
import scipy.sparse as sp
import math
from config_and_print import resolutions, output_directory, software_directory

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

def load_chromosome_sizes(file_path, resolution):
    """
    Load chromosome sizes and calculate matrix dimensions based on resolution.
    
    Parameters:
    - file_path: Path to the chromosome sizes file.
    - resolution: The resolution to divide the chromosome size by.
    
    Returns:
    - chrom_sizes: A dictionary with chromosome names as keys and matrix dimensions as values.
    """
    chrom_sizes = {}
    with open(file_path, 'r') as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = math.ceil(int(size) / resolution)
    return chrom_sizes


def load_sparse_matrix(file_path, shape):
    """
    Load a sparse matrix from a text file.
    
    Parameters:
    - file_path: Path to the sparse matrix file.
    - shape: Tuple indicating the shape of the matrix (rows, columns).
    
    Returns:
    - matrix: A scipy.sparse.csr_matrix representing the matrix.
    """
    rows = []
    cols = []
    data = []
    error_lines = []
    
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            try:
                row, col, value = map(int, line.strip().split())
                if row < 0 or col < 0:
                    raise ValueError(f"Negative row/column index found: row={row}, col={col}")
                rows.append(row)
                cols.append(col)
                data.append(value)
            except ValueError as e:
                error_message = f"Error in file {file_path} at line {i+1}: {str(e)}"
                print(error_message)
                error_lines.append(error_message)

    if error_lines:
        error_log_path = file_path.replace('.txt', '_error_log.txt')
        with open(error_log_path, 'w') as error_log:
            for error in error_lines:
                error_log.write(f"{error}\n")
        print(f"Errors encountered. Details saved in {error_log_path}")

    matrix = sp.csr_matrix((data, (rows, cols)), shape=shape)
    return matrix

def preprocess_sparse_matrix(matrix, regularization=1e-5, threshold=1e-10):
    """
    Preprocess the sparse matrix for better KR normalization.
    """
    matrix.data[matrix.data < threshold] = 0
    matrix.eliminate_zeros()
    
    regularization_matrix = sp.eye(matrix.shape[0], matrix.shape[1]) * regularization
    processed_matrix = matrix + regularization_matrix
    
    return processed_matrix

def save_sparse_matrix(matrix, file_path):
    """
    Save the sparse matrix in the original three-value format (row, col, value).
    
    Parameters:
    - matrix: The scipy.sparse.csr_matrix to save.
    - file_path: The path where the matrix should be saved.
    """
    coo_matrix = matrix.tocoo()
    with open(file_path, 'w') as f:
        for row, col, value in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
            f.write(f"{row}\t{col}\t{value}\n")  

def load_and_process_matrix(prefix, chromosome, shape):
    file_path = f"{output_directory}/hic_{resolution_label}_raw_dir/{chromosome}/{prefix}_{chromosome}.txt"
    
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        print(f"File {file_path} is missing or size zero. Unable to process.")
        return None

    matrix = load_sparse_matrix(file_path, shape)
    processed_matrix = preprocess_sparse_matrix(matrix)
    
    # Save the processed matrix in the original format
    save_sparse_matrix(processed_matrix, file_path.replace('.txt', '_processed.txt'))
    print(f"Processed matrix saved as {file_path.replace('.txt', '_processed.txt')}")
    return processed_matrix

def process_all_chromosomes(chrom_sizes, resolution):
    with open(f"{output_directory}/filtered_bam_list.txt", "r") as file:
        prefixes = file.read().splitlines()

    for chromosome, size in chrom_sizes.items():
        shape = (size, size)  # Matrix shape based on chromosome size
        for prefix in prefixes:
            load_and_process_matrix(prefix, chromosome, shape)

chrom_sizes = load_chromosome_sizes(f"{software_directory}/hg19.autosome.chrom.sizes", resolution)

# Process all chromosomes
process_all_chromosomes(chrom_sizes, resolution)


