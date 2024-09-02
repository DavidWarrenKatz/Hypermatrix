#!/bin/bash

######################################################################################
# This shell script iterates through the BAM directory and creates a filtered list
# of BAM files that have a sufficient number of high-quality reads. Only these files 
# will be used for downstream processing. The parameters for filtering are located 
# in the config file. Symbolic links of all the BAM files that pass filtering are 
# created in the output directory. If the filtered list already exists in the output 
# directory, this step is skipped.
# File: filter_bams.sh
# Author: David Katz (davidkatz02@gmail.com)
######################################################################################

# SAMPLE DATA FOR TEST MODE, WILL UPDATE DEFAULT CONFIGURATON 
# using sample data from https://www.nature.com/articles/s41592-019-0547-z

CONFIG_AND_PRINT_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'config_and_print.py'))")
eval "$(python3 $CONFIG_AND_PRINT_PATH)"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Load samtools module
module load samtools

# Check if the filtered list already exists
if [ ! -f "$filtered_list" ]; then
  # Create or clear the filtered list file
  > "$filtered_list"

  # Count the initial number of BAM files
  initial_bam_count=$(ls "$bam_directory"/sc*.b37.calmd.bam | wc -l)
  echo "Initial number of BAM files: $initial_bam_count"

  # Loop through BAM files and filter based on quality
  for bam_file in "$bam_directory"/sc*.b37.calmd.bam; do
    echo "Evaluating $bam_file"
    
    # Count high-quality reads
    high_quality_reads=$(samtools view -c -q 30 -f 1 -F 1804 "$bam_file")
    
    # Check if the BAM file meets the quality criteria
    if (( high_quality_reads >= min_high_quality_reads )); then
      # Extract the identifier and add it to the filtered list
      identifier=$(basename "$bam_file" .b37.calmd.bam)
      echo "$identifier" >> "$filtered_list"
    fi
  done

  # Count the number of BAM files in the filtered list
  filtered_bam_count=$(wc -l < "$filtered_list")
  echo "Number of BAM files in the filtered list: $filtered_bam_count"
  echo "Filtered list of BAM files has been created."
else
  echo "Filtered list already exists. Skipping filtering step."
fi
