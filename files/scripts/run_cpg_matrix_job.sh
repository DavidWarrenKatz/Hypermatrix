#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=150G
#SBATCH --error=/home/dwk681/workspace/grant/logs/%JgetCpG.err
#SBATCH --output=/home/dwk681/workspace/grant/logs/%JgetCpG.out

# Load the necessary modules for conda
source /etc/profile
source ~/.bashrc

# Activate the conda environment
conda activate multiomics6

python << EOF
import pyBigWig
import numpy as np
from scipy.stats import pearsonr
from scipy import io
import h5py
from scipy.stats import ks_2samp

# Open the bigWig file
path = '/projects/b1198/epifluidlab/david/Genomic_Files/'
data_path="/home/dwk681/workspace/grant/"

bedfile = '/home/dwk681/workspace/Genomic_Files/BedFiles/hg19/ENCFF257GGV_mCpG_hg19_GM12878_WGBS.bed'

# Read the BED file
def read_bed(file_path):
    bed_data = {}
    with open(file_path) as f:
        for line in f:
            if line.startswith("track") or line.startswith("#"):
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            value = float(parts[4]) if len(parts) > 4 else 1.0  # Use the score as the value

            if chrom not in bed_data:
                bed_data[chrom] = []
            bed_data[chrom].append((start, end, value))
    return bed_data

bed_data = read_bed(bedfile)

large_bin_size = 1_000_000  # Define your large bin size
small_bin_size = 10_000  # Define your small bin size

for chrom in bed_data:
    #if chrom not in ('chrX', 'chrY', 'chrM'):
    if chrom in ('chr7', 'chr8', 'chr9', 'chr22'):
        chrom_data = bed_data[chrom]
        chrom_data.sort()
        chrom_length = chrom_data[-1][1]  # Use the end of the last entry as chromosome length
        large_bins = range(0, chrom_length, large_bin_size)
        num_large_bins = len(large_bins)

        integrals = np.zeros((int(large_bin_size / small_bin_size), num_large_bins))

        for i, large_bin_start in enumerate(large_bins):
            for j in range(int(large_bin_size / small_bin_size)):
                small_bin_start = large_bin_start + j * small_bin_size
                small_bin_end = small_bin_start + small_bin_size

                # Ensure that the coordinates are within the chromosome boundaries
                if small_bin_start < 0:
                    small_bin_start = 0
                if small_bin_end > chrom_length:
                    small_bin_end = chrom_length

                # Check if the adjusted small_bin_start is still greater than or equal to small_bin_end
                if small_bin_start >= small_bin_end:
                    continue

                # Get data within the small bin
                small_bin_data = [value for (start, end, value) in chrom_data if start >= small_bin_start and end <= small_bin_end]

                # Calculate the integral of data in the small bin
                if small_bin_data:
                    integral = np.trapz(small_bin_data)
                    integrals[j][i] = integral

        # Create matrices
        num_rows, num_columns = integrals.shape
        p_value_ks_matrix = np.zeros((num_columns, num_columns))
        epsilon = 1e-10  # Define your epsilon value here

        for i in range(num_columns):
            for j in range(i, num_columns):
                if i == j:
                    p_value_ks_matrix[i, j] = 1.0
                else:
                    _, ks_p_value = ks_2samp(integrals[:,i], integrals[:,j])
                    p_value_ks_matrix[i, j] = np.log10(ks_p_value + epsilon)
                    p_value_ks_matrix[j, i] = np.log10(ks_p_value + epsilon)

        # Calculate and save correlation matrices
        correlation_ks_p_value = np.corrcoef(p_value_ks_matrix) + 1

        new_filename = data_path + f'Workspaces/mCpG_{chrom}_res{large_bin_size}_small_{small_bin_size}_ks_p_value_correlation.h5'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=correlation_ks_p_value)

EOF
