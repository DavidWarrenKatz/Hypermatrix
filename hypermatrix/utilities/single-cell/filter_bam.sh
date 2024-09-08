#!/bin/bash

######################################################################################
# This shell script iterates through the BAM directory and creates a filtered list
# of BAM files that have a sufficient number of high-quality reads. Only these files 
# will be used for downstream processing. The parameters for filtering are located 
# in the config file. Symbolic links of all the BAM files that pass filtering are 
# created in the output directory. If the filtered list already exists in the output 
# directory, this step is skipped.
# Additionally, if one of the autosomal chromosomes does not have any high-quality 
# reads, the BAM file will be excluded from the filtered list.
# File: filter_bams.sh
# Author: David Katz (davidkatz02@gmail.com)
######################################################################################

# SAMPLE DATA FOR TEST MODE, WILL UPDATE DEFAULT CONFIGURATON 
# using sample data from https://www.nature.com/articles/s41592-019-0547-z

#K wrote this, I am commenting it out for now, not sure what it is
#CONFIG_AND_PRINT_PATH=$(python3 -c "import pkg_resources; print(pkg_resources.resource_filename('hypermatrix', 'config_and_print.py'))")
#eval "$(python3 $CONFIG_AND_PRINT_PATH)"

eval "$(python3 config_and_print.py)"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Load samtools module
module load samtools

# List of autosomal chromosomes (assuming hg19 or similar)
autosomal_chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

# Check if the filtered list already exists
if [ ! -f "$filtered_list" ]; then
  # Create or clear the filtered list file
  > "$filtered_list"

  # Count the initial number of BAM files
  initial_bam_count=$(ls "$bam_directory"/*.bam | wc -l)
  echo "Initial number of BAM files: $initial_bam_count"

  # Loop through all BAM files and filter based on quality
  for bam_file in "$bam_directory"/*.bam; do
    echo "Evaluating $bam_file"
    
    # Check high-quality reads for each autosomal chromosome
    valid_file=true
    for chr in "${autosomal_chromosomes[@]}"; do
      high_quality_reads_chr=$(samtools view -c -q 30 -f 1 -F 1804 "$bam_file" "$chr")
      if (( high_quality_reads_chr == 0 )); then
        valid_file=false
        echo "No high-quality reads found for $chr in $bam_file. Excluding this file."
        break
      fi
    done

    # Count total high-quality reads in the BAM file
    high_quality_reads_total=$(samtools view -c -q 30 -f 1 -F 1804 "$bam_file")
    
    # Check if the BAM file meets the overall quality criteria and autosomal read check
    if [[ "$valid_file" == true && "$high_quality_reads_total" -ge "$min_high_quality_reads" ]]; then
      # Extract the identifier and add it to the filtered list
      identifier=$(basename "$bam_file" .bam)
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








