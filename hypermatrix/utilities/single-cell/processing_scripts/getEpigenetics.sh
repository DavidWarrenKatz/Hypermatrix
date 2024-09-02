#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=1-1:50:00
#SBATCH --ntasks=6
#SBATCH --mem=100000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetEpigenetic.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetEpigenetic.out

python << EOF
import pyBigWig
import numpy as np
from scipy.stats import pearsonr
from scipy import io
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp


# Open the bigWig file
path = '/projects/b1198/epifluidlab/david/Genomic_Files/'
data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"


# Define bin sizes
large_bin_size = 1000000  
small_bin_size = 1000   

bigwigfile10 = path + 'BigWigFiles/hg19/repressive/ENCFF533NIQ_signal_p-value_hg19_H3K9me3_GM12878.bigWig'   #H3K9me3
#bigwigfile10 = path + 'ENCFF533NIQ_signal_p-value_hg19_H3K9me3_GM12878.bigWig'   #H3K9me3

bw = pyBigWig.open(bigwigfile10)
chromosomes = bw.chroms()

for chrom in chromosomes:
    if chrom not in ('chrX', 'chrY', 'chrM'):
        chrom_length = chromosomes[chrom]
        large_bins = range(0, chrom_length, large_bin_size)
        num_large_bins = len(large_bins)

        integrals = np.zeros((int(large_bin_size / small_bin_size), num_large_bins))

        for i, large_bin_start in enumerate(large_bins):
            for j in range(small_bin_size):
                small_bin_start = large_bin_start + j * small_bin_size
                small_bin_end = small_bin_start + small_bin_size - 1

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
                    integrals[j][i] = integral
                except Exception as e:
                    # Print the exception and relevant information
                    print(f"An error occurred for chromosome {chrom}, large bin {i}, small bin {j}: {e}")

        # create your matrices
        num_rows, num_columns = integrals.shape
        pearson_matrix = np.zeros((num_columns, num_columns)) 
        p_value_pearson_matrix = np.zeros((num_columns, num_columns))  
        ks_matrix = np.zeros((num_columns, num_columns))
        p_value_ks_matrix = np.zeros((num_columns, num_columns))

        epsilon = 1e-10  # Define your epsilon value here

        for i in range(num_columns):
            for j in range(i, num_columns):
                if i == j:
                    p_value_pearson_matrix[i, j] = 1.0  # Diagonal elements should be 1 (correlation with itself is always 1)
                    p_value_ks_matrix[i, j] = 1.0
                    ks_matrix[i,j] = 1.0
                    pearson_matrix[i,j] = 1.0
                else:
                    ks_statistic, ks_p_value = ks_2samp(integrals[:,i], integrals[:,j])
                    
                    if np.var(integrals[:,i]) > 0 and np.var(integrals[:,j]) > 0:
                        pearson_statistic, pearson_p_value = pearsonr(integrals[:,i], integrals[:,j])
                    else:
                        pearson_statistic = 0.0
                        pearson_p_value = 1.0
                    
                    pearson_matrix[i, j] = np.log10(1 + pearson_statistic + epsilon)
                    pearson_matrix[j, i] = np.log10(1 + pearson_statistic + epsilon)
                    p_value_pearson_matrix[i, j] = np.log10(pearson_p_value + epsilon)
                    p_value_pearson_matrix[j, i] = np.log10(pearson_p_value + epsilon)
                    p_value_ks_matrix[i, j] = np.log10(ks_p_value + epsilon)
                    p_value_ks_matrix[j, i] = np.log10(ks_p_value + epsilon)
                    ks_matrix[i, j] = np.log10(ks_statistic + epsilon)
                    ks_matrix[j, i] = np.log10(ks_statistic + epsilon)

                    
        # Save the epigenetic matrix as a .mat file
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_pearsons.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=pearson_matrix)
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_pearsons_p_value.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=p_value_pearson_matrix)            
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_ks.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=ks_matrix)
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_ks_p_value.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=p_value_ks_matrix) 
        
        #take the correlation matrices of all the matrices
        correlation_pearson_statistic = np.corrcoef(pearson_matrix) + 1
        correlation_pearson_p_value = np.corrcoef(p_value_pearson_matrix) + 1
        correlation_ks_statistic = np.corrcoef(ks_matrix) + 1
        correlation_ks_p_value = np.corrcoef(p_value_ks_matrix) + 1

        # Save the correlation epigenetic matrix as a .mat file
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_pearsons_correlation.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=correlation_pearson_statistic)
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_pearsons_p_value_correlation.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=correlation_pearson_p_value)            
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_ks_correlation.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=correlation_ks_statistic)
        new_filename = data_path + f'Workspaces/epigenetic/signal_p_value_H3K9me3_{chrom}_res{large_bin_size}_small_{small_bin_size}_ks_p_value_correlation.mat'
        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('matrix', data=correlation_ks_p_value) 
  

EOF





























