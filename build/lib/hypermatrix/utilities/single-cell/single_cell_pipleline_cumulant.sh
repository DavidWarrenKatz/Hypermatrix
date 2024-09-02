#!/bin/bash
# Description: Run the pipeline for deg three cumulants 
# File: single_cell_pipeline_cumulant.sh
# Plan: The package should be entirely self-contained 

#set -x

echo -e "\n [DEBUG-6]: Entering the Cumulant Script HERE!!!! CUMULANT ARGUMENT SUPPLIED \n"

# Get the correct path to config_and_print.py using Python and pkg_resources
CONFIG_AND_PRINT_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'config_and_print.py'))")

# Check if the path was correctly retrieved
if [ -z "$CONFIG_AND_PRINT_PATH" ]; then
    echo "[ERROR]: Failed to retrieve the path to config_and_print.py"
    exit 1
fi

# Run the Python script and source the output to import the variables
echo "[DEBUG 7]: Using config_and_print.py at $CONFIG_AND_PRINT_PATH"
eval "$(python3 $CONFIG_AND_PRINT_PATH)"
if [ $? -ne 0 ]; then
    echo "[ERROR]: Failed to run config_and_print.py"
    exit 1
fi


# Resolve the path to filter_bam.sh using pkg_resources
FILTER_BAM_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'utilities/single-cell/filter_bam.sh'))")

# Check if the filter_bam.sh script exists and is executable
if [ ! -x "$FILTER_BAM_PATH" ]; then
    echo "[ERROR]: filter_bam.sh does not exist or is not executable at $FILTER_BAM_PATH"
    exit 1
fi

# Execute the filtering script
chmod +x "$FILTER_BAM_PATH"
"$FILTER_BAM_PATH"
if [ $? -ne 0 ]; then
    echo "[ERROR]: filter_bam.sh failed to execute"
    exit 1
fi

# # Execute the filtering script
# chmod +x filter_bam.sh
# ./filter_bam.sh
# if [ $? -ne 0 ]; then
#     echo "[ERROR]: filter_bam.sh failed to execute"
#     exit 1
# fi

# # Generate the interval bed file
# python generate_interval_bed.py

# # Execute the filter bigwig files script
# chmod +x filter_bw.sh
# ./filter_bw.sh

# # Execute the script that filters the filtered hic to contian high quality methy
# #[TO DO: need to update this, right now hard-coded]
# chmod +x update_filted_list.sh
# ./update_filted_list.sh

# # Execute the make pairwise contact script
# chmod +x make_pairwise_contact_format_from_bam.sh
# ./make_pairwise_contact_format_from_bam.sh

# # Execute the preprocessing for hicluster script
# chmod +x preprocess_for_hic_cluster.sh
# ./preprocess_for_hic_cluster.sh

# # Execute the make juicer format script
# # This script forms the raw directory and the 
# #imputed directories
# chmod +x make_juicer_short_format.sh
# ./make_juicer_short_format.sh

# # Make single cells KR normalized if normalization is not NONE
# if [ "$normalization" != "NONE" ]; then
#     echo "Normalization is set to $normalization. Running make_hic_for_normalization.py..."
#     #python make_hic_for_normalization.py
#     chmod +x make_hic_for_normalization.sh
#     ./make_hic_for_normalization.sh
# else
#     echo "Normalization is set to NONE. Skipping normalization step."
# fi

# #Some matrices will be too sparse for KR processing
# #This addiitional step will try to regularize them
# if [ "$normalization" != "NONE" ]; then
#     echo "Normalization is set to $normalization. Running process_sparse_matrices.py ..."
#     python process_sparse_matrices.py
#     chmod +x make_hic_for_normalization_afterprocessing.sh
#     ./make_hic_for_normalization_afterprocessing.sh
# else
#     echo "Normalization is set to NONE. Skipping normalization step."
# fi

# # Source the conda environment setup script
# source /software/miniconda3/4.12.0/etc/profile.d/conda.sh
# # Define environment names
# hypermatrix_env=hypermatrix
# # Activate the desired conda environment
# conda activate $hypermatrix_env

# # Make the methylation matrices
# python make_methy_matrices.py

# # Make the hic matrices
# python make_hic_matrices.py

# # Make the hic cumulants if cumulant is set to True
# #[TO DO: make this file]
# if [ "$impute" = "True" ]; then
# python make_hic_matrices_imputed.py
# fi

# python make_hic_cumulant.py

# # Make the hic cumulants if cumulant is set to True
# #[TO DO: make this file]
# if [ "$cumulant" = "True" ]; then
# python make_methy_cumulant.py
# fi

# #Make the combined tensor for each cell
# #[TO DO: need to decide if imputation and emphasis and merging will happen]i
# #Right now, emphasis, correlation, shift
# python make_combined_cumulant_tensor.py



# comment

# #Check to make sure tensorlab is available
# # Define the paths
# zip_file="../../bin/softwarefiles/tensorlab/tensorlab_2016-03-28.zip"
# echo "DEBUG-7 expecting the following zipfile $zip_file"
# unzip_dir="../../bin/softwarefiles/tensorlab/"

# # Check if the zip file exists
# if [ ! -f "$zip_file" ]; then
#   echo "Zip file does not exist."
#   exit 1
# fi

# # Get the list of files in the zip archive
# zip_files=$(unzip -l "$zip_file" | awk 'NR>3 {print $4}' | sed '$d')

# # Function to check if all files from the zip are already in the directory
# all_files_exist() {
#   for file in $zip_files; do
#     if [ ! -f "$unzip_dir/$file" ]; then
#       return 1
#     fi
#   done
#   return 0
# }

# # Check if all files exist
# if all_files_exist; then
#   echo "All files from the zip are already present in the directory. Skipping unzip."
# else
#   echo "Unzipping the file..."
#   unzip "$zip_file" -d "$unzip_dir"
#   echo "Unzipping completed."
# fi

# #compute the AB compartments for each cell, display results
# # Load the necessary modules
# module load matlab/r2022b

# # Execute the MATLAB script
# matlab -nodisplay -r "run('get_AB_single_cell_structured_data_cumulant.m'); exit;"

# # Download the dark regions file if it doesn't already exist
# dark_regions_file="../../bin/softwarefiles/dark_regions_hg19.bigWig"
# if [ ! -f "$dark_regions_file" ]; then
#     wget $dark_regions_hg19_url -O "$dark_regions_file"
# fi

# #Create dark bins file if it does not already exist
# python create_dark_bins.py

# chmod +x get_eigenvectors_bulk.sh 
# ./get_eigenvectors_bulk.sh

# #Execute script for making AB calls
# python make_AB_compartments.py

# #Make image of AB calls
# #This part is specfic to IMR90 GM12878 experiment
# python make_AB_compartment_image.py

# python make_all_cells_tensor.py

# # Load the necessary modules
# module load matlab/r2022b
# # Execute the MATLAB script
# matlab -nodisplay -r "run('get_cell_type_factors_all_cell.m'); exit;"

# python make_cell_type_factors.py

# module load matlab/r2022b
# matlab -nodisplay -r "run('cell_type_factors_seperate_modalities.m'); exit;"

# #Execute script to find differential bins
# #python find_differential_bins.py

# #analyze AB compartment calls

# #cluster cells

# #find highly variable AB compartment region from bulk data

# #display clustered cells colored by their AB compartment 
# #value for highly variable region

# #if cumulnat flag is TRUE, need to repeat the clustering and calls 
# #with 3rd degree cumulants

