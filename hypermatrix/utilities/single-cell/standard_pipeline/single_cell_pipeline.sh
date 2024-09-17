#!/bin/bash
# File: single_cell_pipeline.sh

##########################################################################
# Description: The following script converts bam files into hic matrices #
# and methylation matrices organized in folders for each chromosome.     #
# This pipeline also computes comparts and clustering and outputs        #
# images for comparison.                                                 #
##########################################################################



# 1. load export_config.py 
EXPORT_CONFIG_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'export_config.py'))")
MAKE_METHY_MATRICES_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'make_methy_matrices.py'))")
MAKE_HIC_MATRICES_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'make_hic_matrices.py'))")
MAKE_COMBINED_METHY_HIC_MATRICES_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'make_combined_methy_hic_tensor_single_cell.py'))")



# 1. Load scripts from package 
# Check if the path was correctly retrieved
# Write function to resolve script tools loop through 
if [ -z "$EXPORT_CONFIG_PATH" ]; then
    echo "[ERROR]: Failed to retrieve the path to export_config.py"
    exit 1
else
    echo "[LOG]: successfully resolved path for export_config.py \n"
    echo "\n PATH: ${EXPORT_CONFIG_PATH} s\n"
fi


# Check if the path was correctly retrieved
if [ -z "$MAKE_METHY_MATRICES_PATH" ]; then
    echo "[ERROR]: Failed to retrieve the path to export_config.py"
    exit 1
else
    echo "\n [LOG]: successfully resolved path for export_config.py \n"
    echo "\n PATH: ${MAKE_METHY_MATRICES_PATH} s\n"
fi

# Check if the path was correctly retrieved
if [ -z "$MAKE_HIC_MATRICES_PATH" ]; then
    echo "[ERROR]: Failed to retrieve the path to export_config.py"
    exit 1
else
    echo "[LOG]: successfully resolved path for export_config.py \n"
    echo "\n PATH: ${MAKE_HIC_MATRICES_PATH} s\n"
fi


# Check if the path was correctly retrieved
if [ -z "$MAKE_COMBINED_METHY_HIC_MATRICES_PATH" ]; then
    echo "[ERROR]: Failed to retrieve the path to export_config.py"
    exit 1
else
    echo "[LOG]: successfully resolved path for export_config.py \n"
    echo "\n PATH: ${MAKE_COMBINED_METHY_HIC_MATRICES_PATH} s\n"
fi

# Import parameters from export_config.py 
# Load as variables into work space
eval "$(python3 "${EXPORT_CONFIG_PATH}")"



# # Get the correct path to config_and_print.py using Python and pkg_resources
# CONFIG_AND_PRINT_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'config_and_print.py'))")




# # Get the directory where this script (single_cell_pipeline.sh) is located
# SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# echo "[LOG] -- Script directory: $SCRIPT_DIR"

# # Import the parameters from config.py (relative to the script's directory)
# eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

# # Make the methylation matrices
# # FIX with AB2
# python $SCRIPT_DIR/make_methy_matrices.py


# # Fix with AB3
# # Make the hic matrices
# python $SCRIPT_DIR/make_hic_matrices.py

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
