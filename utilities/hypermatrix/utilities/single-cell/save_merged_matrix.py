import os
import glob
import numpy as np
import argparse
import re

# Parse command-line arguments for resolution, label, and base directory
parser = argparse.ArgumentParser(description='Save merged matrix for different resolutions.')
parser.add_argument('--label', type=str, help='Resolution label, e.g., "1Mb"')
parser.add_argument('--base_dir', type=str, help='Base directory for input and output files')
args = parser.parse_args()

label = args.label
base_dir = os.path.join(args.base_dir, f'hicluster_{label}_impute_dir/merged/')

# Pattern to match all .npy files in the merged directory for the specified resolution
pattern = os.path.join(base_dir, f'pad1_std1_rp0.5_sqrtvc_chr*.cpgcomp.npy')

# Find all .npy files matching the pattern
npy_files = glob.glob(pattern)

# Function to extract the chromosome number for sorting
def extract_chr_num(filename):
    match = re.search(r'chr(\d+|\X)', filename)
    if match:
        chr_num = match.group(1)
        return int(chr_num) if chr_num.isdigit() else 25 if chr_num == 'X' else 0
    return 0  # Default case for sorting

# Sort files based on the chromosome number
npy_files_sorted = sorted(npy_files, key=extract_chr_num)

# Load each .npy file, drop the last bin, and append to the list
loaded_arrays = [np.delete(np.load(file), -1, axis=1) for file in npy_files_sorted]

# Concatenate all arrays along the first axis (rows) to match desired dimensions
large_matrix = np.concatenate(loaded_arrays, axis=1)
large_matrix = large_matrix.T

# Verify the shape of the large_matrix
print(f"The shape of the large matrix is: {large_matrix.shape}")

# Specify the output file path for saving the matrix, adjusted for the specified resolution
output_dir = os.path.join(args.base_dir, 'matrices', label)
os.makedirs(output_dir, exist_ok=True)
output_file_path = os.path.join(output_dir, f'b37.autosome.{label}_interval.add_value.hic.GM_IMR90_153_samples.npy')

# Save the matrix to a file
np.save(output_file_path, large_matrix)
