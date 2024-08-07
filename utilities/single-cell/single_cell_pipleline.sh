#!/bin/bash
# Description: The following script converts bam files into hic matrices 
# and methylation matrices organized in folders for each chromosome.
# This pipeline also computes comparts and clustering and outputs
# images for comparison.
# [TO DO: everything runs in sequence, rewrite to make it run in parallel]

# Run the Python script and source the output to import the variables
# each individual file still needs to import these variables as well
eval "$(python3 config_and_print.py)"

# Execute the filtering script
chmod +x filter_bam.sh
./filter_bam.sh

# Execute the filter bigwig files script
chmod +x filter_bw.sh
./filter_bw.sh

# Execute the script that filters the filtered hic to contian high quality methy
#[TO DO: need to update this, right now hard-coded]
chmod +x update_filted_list.sh
./update_filted_list.sh

# Execute the make pairwise contact script
chmod +x make_pairwise_contact_format_from_bam.sh
./make_pairwise_contact_format_from_bam.sh

# Execute the preprocessing for hicluster script
chmod +x preprocess_for_hic_cluster.sh
./preprocess_for_hic_cluster.sh

# Execute the make juicer format script
# This script forms the raw directory and the 
#imputed directories
chmod +x make_juicer_short_format.sh
./make_juicer_short_format.sh

# Execute the compartment calling scripts from scHICluster
# [TO DO: need to get this working]
if [ "$cluster_compartments" = "True" ]; then
chmod +x compartment_calling.sh
./compartment_calling.sh
fi

# Source the conda environment setup script
source /software/miniconda3/4.12.0/etc/profile.d/conda.sh
# Define environment names
hypermatrix_env=hypermatrix
# Activate the desired conda environment
conda activate $hypermatrix_env

# Make the methylation matrices
python make_methy_matrices.py

# Make the hic matrices
python make_hic_matrices.py

# Make the hic cumulants if cumulant is set to True
#[TO DO: make this file]
if [ "$cluster_compartments" = "True" ]; then
python make_hic_cumulant.py
fi

#Make the combined tensor for each cell
#[TO DO: need to decide if imputation and emphasis and merging will happen]
python make_combined_methy_hic_tensor_single_cell.py


#compute the AB compartments for each cell, display results

#cluster cells

#find highly variable AB compartment region from bulk data

#display clustered cells colored by their AB compartment 
#value for highly variable region

#if cumulnat flag is TRUE, need to repeat the clustering and calls 
#with 3rd degree cumulants
