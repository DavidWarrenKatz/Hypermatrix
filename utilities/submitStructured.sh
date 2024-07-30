#!/bin/bash

# Activate the conda environment
source /etc/profile
source ~/.bashrc
conda activate hypermatrix

# Load the necessary modules
module load matlab/r2022b

# Execute the MATLAB script
matlab -nodisplay -r "run('getStructuredData.m'); exit;"
