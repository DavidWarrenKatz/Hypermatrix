import matplotlib.pyplot as plt
import h5py
import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, NMF
from sklearn.manifold import TSNE
import tensorly as tl
from tensorly.decomposition import non_negative_parafac
from config_and_print import resolutions, output_directory, filtered_list, normalization

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

# Other utilities
base_directory = output_directory

def normalize_matrix_columns(A):
    """
    Normalize each column of the matrix so that each has a norm of one.
    
    Parameters:
    - A: a 2D NumPy array (matrix)
    
    Returns:
    - normalized_A: a matrix where each column of A has been divided by its L2 norm
    """
    # Calculate the L2 norm for each column
    column_norms = np.linalg.norm(A, axis=0)
    
    # Avoid division by zero
    if np.any(column_norms == 0):
        raise ValueError("One or more columns have zero norm. Cannot normalize those columns.")
    
    # Normalize each column by its norm
    normalized_A = A / column_norms
    
    return normalized_A

def calculate_matrix(file_path):
    # Initialize an empty list to store the data
    data = []

    # Read the compressed file using gzip
    with gzip.open(file_path, 'rt') as file:
        for line in file:
            # Split the line into fields
            fields = line.strip().split('\t')

            # Extract relevant entries for calculations
            numerator_indices = range(6, len(fields), 2)
            denominator_indices = range(7, len(fields), 2)

            row_data = []
            for num_idx, denom_idx in zip(numerator_indices, denominator_indices):
                numerator = float(fields[num_idx])
                denominator = float(fields[denom_idx])

                # Check if denominator is zero, set entry to 0, else perform the division
                row_data.append(0 if denominator == 0 else numerator / denominator)

            data.append(row_data)

    # Convert the list of lists into a NumPy array
    data_matrix = np.array(data)

    return data_matrix

file_path = base_directory + f'b37.autosome.{resolution_label}_interval.add_value.methy.bed.gz'
methylation_matrix = normalize_matrix_columns(calculate_matrix(file_path))
print(methylation_matrix.shape)

# Define the path to the file
file_path = f'{output_directory}/filtered_bam_list.txt'

# Initialize an empty list to store the prefixes
prefixes = []

try:
    # Open and read the file with a different encoding
    with open(file_path, 'r', encoding='ISO-8859-1') as file:
        for line in file:
            # Strip the newline character and append to the prefixes list
            prefixes.append(line.strip())
except UnicodeDecodeError:
    print("Failed to decode file with ISO-8859-1. Trying another method.")

# Print the list of prefixes
print(len(prefixes))

file_list = prefixes

chromosomes_info = {
    '1': 249250621,
    '2': 243199373,
    '3': 198022430,
    '4': 191154276,
    '5': 180915260,
    '6': 171115067,
    '7': 159138663,
    '8': 146364022,
    '9': 141213431,
    '10': 135534747,
    '11': 135006516,
    '12': 133851895,
    '13': 115169878,
    '14': 107349540,
    '15': 102531392,
    '16': 90354753,
    '17': 81195210,
    '18': 78077248,
    '19': 59128983,
    '20': 63025520,
    '21': 48129895,
    '22': 51304566,
}

# Set the desired rank
rank = 2

# Perform Non-negative Matrix Factorization (NMF)
model = NMF(n_components=rank, init='random', random_state=0)
W = model.fit_transform(methylation_matrix_1Mb)
H = model.components_

# Check the shapes of W and H to confirm the decomposition
print("Shape of W:", W.shape)
print("Shape of H:", H.shape)


def get_bins_per_chromosome(chromosomes_info, resolution):
    """
    Calculate the number of bins for each chromosome based on the resolution.
    
    Parameters:
    - chromosomes_info: Dictionary with chromosome names as keys and their lengths as values.
    - resolution: The resolution of the bins.
    
    Returns:
    - bins_per_chromosome: Dictionary with chromosome names as keys and number of bins as values.
    """
    bins_per_chromosome = {}
    for chrom, size in chromosomes_info.items():
        bins = size // resolution
        bins_per_chromosome[chrom] = bins
    return bins_per_chromosome

def split_genomic_factors(W, bins_per_chromosome):
    """
    Split the genomic factors matrix W into individual chromosome vectors.
    
    Parameters:
    - W: Genomic factors matrix (shape: total_bins x rank)
    - bins_per_chromosome: Dictionary with chromosome names as keys and number of bins as values.
    
    Returns:
    - chromosome_vectors: Dictionary with chromosome names as keys and their corresponding factor vectors as values.
    """
    chromosome_vectors = {}
    start_idx = 0
    
    for chrom, bins in bins_per_chromosome.items():
        end_idx = start_idx + bins
        chromosome_vectors[chrom] = W[start_idx:end_idx, :]
        start_idx = end_idx
    
    return chromosome_vectors

# Calculate the number of bins per chromosome
bins_per_chromosome = get_bins_per_chromosome(chromosomes_info, resolution)

# Split the genomic factors into individual chromosome vectors
chromosome_vectors = split_genomic_factors(W, bins_per_chromosome)

# Print the shape of each chromosome vector to confirm
for chrom, vector in chromosome_vectors.items():
    print(f"Chromosome {chrom}: {vector.shape}")



import os

# Define the output directory
output_dir = os.path.join(base_directory, "chromosome_vectors")
os.makedirs(output_dir, exist_ok=True)

# Function to save each vector to a file
def save_vectors_to_files(chromosome_vectors, output_dir):
    for chrom, vector in chromosome_vectors.items():
        file_path = os.path.join(output_dir, f"chromosome_{chrom}_vectors.txt")
        np.savetxt(file_path, vector, fmt='%.6f')
        print(f"Saved {file_path}")

# Function to save the genome-wide vectors
def save_genome_wide_vectors(W, output_dir):
    genome_wide_file = os.path.join(output_dir, "genome_wide_vectors.txt")
    np.savetxt(genome_wide_file, W, fmt='%.6f')
    print(f"Saved genome-wide vectors to {genome_wide_file}")

# Save chromosome vectors
save_vectors_to_files(chromosome_vectors, output_dir)

# Save genome-wide vectors
save_genome_wide_vectors(W, output_dir)

