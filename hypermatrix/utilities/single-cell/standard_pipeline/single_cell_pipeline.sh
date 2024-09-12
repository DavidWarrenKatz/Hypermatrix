#!/bin/bash
# Description: The following script converts bam files into hic matrices 
# and methylation matrices organized in folders for each chromosome.
# This pipeline also computes comparts and clustering and outputs
# images for comparison.
# [TO DO: everything runs in sequence, rewrite to make it run in parallel]
# [TO DO: need to replicate the KR normalization from bulk]
# [TO DO: add sex chromosomes and contigs and mitochondrial DNA analysis]
# sex chromosomes will be important for differentiat chromosomes command.
# Right now, o/e is only done in bulk data, no single cell pipeline
# re-write so that the code does not use hicluster to make bin 0-indexed short format from .hic files
# the matrices produced by scHIcluster do not seem quite right. They are missing the diagonal for example

# Get the directory where this script (single_cell_pipeline.sh) is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Import the parameters from config.py (relative to the script's directory)
eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

# Execute the filtering script (relative to the script's directory)
chmod +x $SCRIPT_DIR/filter_bam.sh
$SCRIPT_DIR/filter_bam.sh

# Generate the interval bed file
python $SCRIPT_DIR/generate_interval_bed.py

# Execute the filter bigwig files script
chmod +x $SCRIPT_DIR/filter_bw.sh
$SCRIPT_DIR/filter_bw.sh

# Execute the script that filters the filtered hic to contian high quality methy
#[TO DO: need to update this, right now hard-coded]
chmod +x $SCRIPT_DIR/update_filted_list.sh
$SCRIPT_DIR/update_filted_list.sh

# Execute the make pairwise contact script
chmod +x $SCRIPT_DIR/make_pairwise_contact_format_from_bam.sh
$SCRIPT_DIR/make_pairwise_contact_format_from_bam.sh

# Execute the preprocessing for hicluster script
chmod +x $SCRIPT_DIR/preprocess_for_hic_cluster.sh
$SCRIPT_DIR/preprocess_for_hic_cluster.sh

# Execute the make juicer format script
# This script forms the raw directory and the 
#imputed directories
chmod +x $SCRIPT_DIR/make_juicer_short_format.sh
$SCRIPT_DIR/make_juicer_short_format.sh

# Make the methylation matrices
python $SCRIPT_DIR/make_methy_matrices.py

# Make the hic matrices
python $SCRIPT_DIR/make_hic_matrices.py

#Make the combined tensor for each cell
#[TO DO: need to decide if imputation and emphasis and merging will happen]
#Right now, emphasis, correlation, shift
python $SCRIPT_DIR/make_combined_methy_hic_tensor_single_cell.py

#Check to make sure tensorlab is available
# Define the paths
zip_file="$SCRIPT_DIR/../../../../src/softwarefiles/tensorlab/tensorlab_2016-03-28.zip"
unzip_dir="$SCRIPT_DIR/../../../../src/softwarefiles/tensorlab/"

# Check if the zip file exists
if [ ! -f "$zip_file" ]; then
  echo "Zip file does not exist."
  exit 1
fi

# Get the list of files in the zip archive
zip_files=$(unzip -l "$zip_file" | awk 'NR>3 {print $4}' | sed '$d')

# Function to check if all files from the zip are already in the directory
all_files_exist() {
  for file in $zip_files; do
    if [ ! -f "$unzip_dir/$file" ]; then
      return 1
    fi
  done
  return 0
}

# Check if all files exist
if all_files_exist; then
  echo "All files from the zip are already present in the directory. Skipping unzip."
else
  echo "Unzipping the file..."
  unzip "$zip_file" -d "$unzip_dir"
  echo "Unzipping completed."
fi

#compute the AB compartments for each cell, display results
# Load the necessary modules
module load matlab/r2022b

# Execute the MATLAB script
matlab -nodisplay -r "run('$SCRIPT_DIR/get_AB_single_cell_structured_data.m'); exit;"

# Download the dark regions file if it doesn't already exist
dark_regions_file="$SCRIPT_DIR/../../../src/softwarefiles/dark_regions_$reference_genome.bigWig"
if [ ! -f "$dark_regions_file" ]; then
    wget $dark_regions_hg_url -O "$dark_regions_file"
fi

#Create dark bins file if it does not already exist
python $SCRIPT_DIR/create_dark_bins.py

chmod +x $SCRIPT_DIR/get_eigenvectors_bulk.sh 
$SCRIPT_DIR/get_eigenvectors_bulk.sh

#Execute script for making AB calls
python $SCRIPT_DIR/make_AB_compartments.py

#Make image of AB calls
#This part is specfic to IMR90 GM12878 experiment
python $SCRIPT_DIR/make_AB_compartment_image.py

<<comment
python make_all_cells_tensor.py

# Load the necessary modules
module load matlab/r2022b
# Execute the MATLAB script
matlab -nodisplay -r "run('get_cell_type_factors_all_cell.m'); exit;"

python make_cell_type_factors.py

module load matlab/r2022b
matlab -nodisplay -r "run('cell_type_factors_seperate_modalities.m'); exit;"

#Execute script to find differential bins
#python find_differential_bins.py

#analyze AB compartment calls

#cluster cells

#find highly variable AB compartment region from bulk data

#display clustered cells colored by their AB compartment 
#value for highly variable region

comment
