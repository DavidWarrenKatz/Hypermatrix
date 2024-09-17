#!/bin/bash
# File: single_cell_pipeline.sh

##########################################################################
# Description: The following script converts BAM files into Hi-C matrices #
# and methylation matrices organized in folders for each chromosome.      #
# This pipeline also computes compartments and clustering and outputs     #
# images for comparison.                                                  #
##########################################################################

# Enable robust debugging
set -x  # Prints commands and their arguments before execution

# Function to check if a file path is resolved correctly
check_path() {
    local script_name=$1
    local path=$2

    echo "[DEBUG]: Checking if the path for $script_name is valid..."
    if [ ! -f "$path" ]; then
        echo "[ERROR]: Failed to retrieve the path to $script_name (Path: $path)"
        exit 1
    else
        echo "[LOG]: Successfully resolved path for $script_name"
        echo " PATH: $path"
    fi
}

# Function to execute a Python script and check its exit status
run_script() {
    local script_path=$1
    local script_desc=$2

    echo "[LOG]: Running $script_desc (Path: $script_path)"
    python "$script_path"
    
    local exit_code=$?
    echo "[DEBUG]: Exit code for $script_desc is $exit_code"
    
    if [ $exit_code -eq 0 ]; then
        echo "[LOG]: $script_desc completed successfully"
    else
        echo "[ERROR]: $script_desc failed to execute with exit code $exit_code"
        exit 1
    fi
}

# Debug: Start of script
echo "[DEBUG]: Starting single-cell pipeline script"

# 1. Resolve script paths by constructing them relative to the hypermatrix package
echo "[DEBUG]: Resolving script paths using known directory structure"

# Get the directory of the hypermatrix package
HYPERMATRIX_DIR=$(python3 -c "
import os
import hypermatrix
print(os.path.dirname(hypermatrix.__file__))
")

# Construct paths to the scripts
EXPORT_CONFIG_PATH="$HYPERMATRIX_DIR/export_config.py"
MAKE_METHY_MATRICES_PATH="$HYPERMATRIX_DIR/utilities/single-cell/standard_pipeline/make_methy_matrices.py"
MAKE_HIC_MATRICES_PATH="$HYPERMATRIX_DIR/utilities/single-cell/standard_pipeline/make_hic_matrices.py"
MAKE_COMBINED_METHY_HIC_MATRICES_PATH="$HYPERMATRIX_DIR/utilities/single-cell/standard_pipeline/make_combined_methy_hic_tensor_single_cell.py"

# Debug: Print resolved paths
echo "[DEBUG]: Resolved paths:"
echo "EXPORT_CONFIG_PATH: $EXPORT_CONFIG_PATH"
echo "MAKE_METHY_MATRICES_PATH: $MAKE_METHY_MATRICES_PATH"
echo "MAKE_HIC_MATRICES_PATH: $MAKE_HIC_MATRICES_PATH"
echo "MAKE_COMBINED_METHY_HIC_MATRICES_PATH: $MAKE_COMBINED_METHY_HIC_MATRICES_PATH"

# Check paths using the check_path function
check_path "export_config.py" "$EXPORT_CONFIG_PATH"
check_path "make_methy_matrices.py" "$MAKE_METHY_MATRICES_PATH"
check_path "make_hic_matrices.py" "$MAKE_HIC_MATRICES_PATH"
check_path "make_combined_methy_hic_tensor_single_cell.py" "$MAKE_COMBINED_METHY_HIC_MATRICES_PATH"

# 3. Run the scripts and check exit status
run_script "$EXPORT_CONFIG_PATH" "Export configuration"
run_script "$MAKE_METHY_MATRICES_PATH" "Make methylation matrices"
run_script "$MAKE_HIC_MATRICES_PATH" "Make Hi-C matrices"
run_script "$MAKE_COMBINED_METHY_HIC_MATRICES_PATH" "Make combined methylation and Hi-C tensor"




# echo "[LOG]: preparing to export configuration"
# eval "$(python3 "${EXPORT_CONFIG_PATH}")" 

# # -eq 0 || echo "[ERROR]: failed to retrieve output configurations"

# echo "[LOG]: Generating methylation matrix generation"
# python $MAKE_METHY_MATRICES_PATH

# # Ref: http://homer.ucsd.edu/homer/interactions/HiCmatrices.html
# echo "[LOG]: Generating HIC matrices methylation matrix generation"
# python $MAKE_METHY_MATRICES_PATH

# echo "[LOG]: Combining Methyl + HiC Single Cell in Tensor"
# python $MAKE_METHY_MATRICES_PATH



# #Make the combined tensor for each cell
# #[TO DO: need to decide if imputation and emphasis and merging will happen]
# #Right now, emphasis, correlation, shift
# python $SCRIPT_DIR/make_combined_methy_hic_tensor_single_cell.py

# #Check to make sure tensorlab is available
# # Define the paths
# zip_file="$SCRIPT_DIR/../../../../src/softwarefiles/tensorlab/tensorlab_2016-03-28.zip"
# unzip_dir="$SCRIPT_DIR/../../../../src/softwarefiles/tensorlab/"

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
# matlab -nodisplay -r "run('$SCRIPT_DIR/get_AB_single_cell_structured_data.m'); exit;"

# # Download the dark regions file if it doesn't already exist
# dark_regions_file="$SCRIPT_DIR/../../../src/softwarefiles/dark_regions_$reference_genome.bigWig"
# if [ ! -f "$dark_regions_file" ]; then
#     wget $dark_regions_hg_url -O "$dark_regions_file"
# fi

# #Create dark bins file if it does not already exist
# python $SCRIPT_DIR/create_dark_bins.py

# chmod +x $SCRIPT_DIR/get_eigenvectors_bulk.sh 
# $SCRIPT_DIR/get_eigenvectors_bulk.sh

# #Execute script for making AB calls
# python $SCRIPT_DIR/make_AB_compartments.py

# #Make image of AB calls
# #This part is specfic to IMR90 GM12878 experiment
# python $SCRIPT_DIR/make_AB_compartment_image.py

# <<comment
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

# comment




# # [TO DO: everything runs in sequence, rewrite to make it run in parallel]
# # [TO DO: need to replicate the KR normalization from bulk]
# # [TO DO: add sex chromosomes and contigs and mitochondrial DNA analysis]
# # sex chromosomes will be important for differentiat chromosomes command.
# # Right now, o/e is only done in bulk data, no single cell pipeline
# # re-write so that the code does not use hicluster to make bin 0-indexed short format from .hic files
# # the matrices produced by scHIcluster do not seem quite right. They are missing the diagonal for example

# # alternative method is to resolve path using pkg_resources
# # all packages should be resolved relative to python3.9/site-packages/hypermatrix 
