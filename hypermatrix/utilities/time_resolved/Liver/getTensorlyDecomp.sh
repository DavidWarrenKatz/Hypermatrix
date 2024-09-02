#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=3-9:00:00
#SBATCH --ntasks=5
#SBATCH --mem=150000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetStructuredDataLog.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetStructuredDataLog.out


# Set the path for data storage
data_path="/home/dwk681/workspace/CRA004660/Liver/"


python << EOF
path = '$data_path';

import tensorly as tl
import numpy as np
from tensorly.decomposition import non_negative_parafac
import matplotlib.pyplot as plt
import h5py

# Step 1: Load the HDF5 file
with h5py.File('/home/dwk681/workspace/CRA004660/Liver/tensor_normalized_filtered_339_500_32285.h5', 'r') as hf:
    # Step 2: Load the tensor
    tensor = hf['tensor'][:]

    # Step 3: Load the labels
    labels = hf['labels'][:]

# Step 4: Convert labels back to a list of strings
labels = labels.astype(str).tolist()


# Perform non-negative PARAFAC decomposition
rank = 4  # Specify the rank of the decomposition
factors = non_negative_parafac(tensor, rank=rank, n_iter_max=300, tol=1e-8, init='random', verbose=0)

with h5py.File(f'/home/dwk681/workspace/CRA004660/Liver/factors_{rank}_300iteration_normalized_filtered_{tensor.shape[0]}_{tensor.shape[1]}_{tensor.shape[2]}.h5', 'w') as hf:
    hf.create_dataset(f'factors', data=factors)
    hf.create_dataset('labels', data=np.array(labels, dtype='S'))

print("Factors saved as an HDF5 file")

EOF
