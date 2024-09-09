#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=1:48:00
#SBATCH --ntasks=5
#SBATCH --mem=100000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/IMR90/logs/%JgetOtherData.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/IMR90/logs/%JgetOtherData.out

# Load the necessary modules for conda
source /etc/profile
source ~/.bashrc

# Activate the conda environment
conda activate multiomics6

# Set the path for data storage
data_path="/projects/b1198/epifluidlab/david/GSE63525/IMR90/"

python << FILE
import numpy as np
import scipy.io
import h5py
path = '$data_path'
resolutions = [1000000]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
data_types = ['observed']

def rows_with_ninety_percent_zeros(matrix):
    rows_removed_indices = []
    num_cols = len(matrix[0])
    threshold = 0.95 * num_cols  # Threshold for 90% zeros in a row

    for idx, row in enumerate(matrix):
        zero_count = sum(1 for value in row if value == 0)
        if zero_count >= threshold:
            rows_removed_indices.append(idx)

    return rows_removed_indices

for resolution in resolutions:
    for chromosome in chromosomes:
        for data in data_types:
            mat_data = path + f'Workspaces/individual/ch{chromosome}_res{resolution}_{data}_NONE.mat'
            file = h5py.File(mat_data, 'r')
            matrix = file['matrix'][:]

            # Get the indices of removed rows
            removed_indices = rows_with_ninety_percent_zeros(matrix)
            # Convert the indices to a NumPy array
            result_array = np.array(removed_indices)
            # Save the indices to a MATLAB .mat file
            scipy.io.savemat(path + f'Workspaces/individual/ch{chromosome}_res{resolution}_{data}_NONE_removedRows.mat', {'removed_rows_indices': result_array})
FILE
