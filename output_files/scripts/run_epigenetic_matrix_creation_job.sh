#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=150G
#SBATCH --error=/home/dwk681/workspace/grant/logs/%JgetEpigenetic.err
#SBATCH --output=/home/dwk681/workspace/grant/logs/%JgetEpigenetic.out

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

bigwigfile = path + 'BigWigFiles/hg19/active/ENCFF372EXQ_hg19_H3K4me1_GM12878.bigWig'   
bw = pyBigWig.open(bigwigfile)
chromosomes = bw.chroms()

large_bin_size = 1_000_000
small_bin_size = 10_000

for chrom in chromosomes:
    if chrom not in ('chrX', 'chrY', 'chrM'):
        chrom_length = chromosomes[chrom]
        large_bins = range(0, chrom_length, large_bin_size)
        num_large_bins = len(large_bins)

        integrals = np.zeros((large_bin_size // small_bin_size, num_large_bins))

        for i, large_bin_start in enumerate(large_bins):
            for j in range(large_bin_size // small_bin_size):
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
                    
                try:
                    # Get data within the small bin
                    small_bin_data = bw.values(chrom, int(small_bin_start), int(small_bin_end))

                    # Calculate the integral of data in the small bin
                    integral = np.trapz(small_bin_data)

                    # Store the integral in the corresponding list
                    integrals[j, i] = integral
                except Exception as e:
                    # Print the exception and relevant information
                    print(f"An error occurred for chromosome {chrom}, large bin {i}, small bin {j}: {e}")

        # create your matrices
        num_rows, num_columns = integrals.shape
        p_value_ks_matrix = np.zeros((num_columns, num_columns))
        epsilon = 1e-10  # Define your epsilon value here
        
        for i in range(num_columns):
            for j in range(i, num_columns):
                if i == j:
                    p_value_ks_matrix[i, j] = 1.0
                else:
                    _, ks_p_value = ks_2samp(integrals[:, i], integrals[:, j])

                    p_value_ks_matrix[i, j] = np.log10(ks_p_value + epsilon)
                    p_value_ks_matrix[j, i] = np.log10(ks_p_value + epsilon)
                    
        #take the correlation matrices of all the matrices
        correlation_ks_p_value = np.corrcoef(p_value_ks_matrix) + 1
        
        # Save the correlation epigenetic matrix as a .mat file
        new_filename = data_path + f'Workspaces/signal_p_value_H3K4me1_{chrom}_res{large_bin_size}_small_{small_bin_size}_ks_p_value_correlation.h5'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=correlation_ks_p_value)

EOF
