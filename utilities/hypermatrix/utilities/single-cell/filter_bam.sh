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

# Load configuration and print variables
eval "$(python3 config_and_print.py)"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Load samtools module
module load samtools

# Function to retrieve autosomal chromosomes from the BAM file
get_autosomal_chromosomes() {
  local bam_file="$1"
  # Extract chromosome names from BAM file using idxstats
  chromosomes=$(samtools idxstats "$bam_file" | cut -f 1)
  
  # Filter autosomal chromosomes (excluding sex chromosomes and mitochondrial DNA)
  autosomal_chromosomes=()
  for chr in $chromosomes; do
    # Only include chromosomes that look like autosomes (i.e., numbers or chr<number>)
    if [[ "$chr" =~ ^(chr)?[0-9]+$ ]]; then
      autosomal_chromosomes+=("$chr")
    fi
  done
  
  echo "${autosomal_chromosomes[@]}"
}

# Check if the filtered list already exists
if [ ! -f "$filtered_list" ]; then
  # Create a temporary filtered list file
  temp_filtered_list=$(mktemp)
  
  # Count the initial number of BAM files
  initial_bam_count=$(ls "$bam_directory"/*.bam | wc -l)
  echo "Initial number of BAM files: $initial_bam_count"

  # Loop through all BAM files and filter based on quality
  for bam_file in "$bam_directory"/*.bam; do
    echo "Evaluating $bam_file"

    # Retrieve the list of autosomal chromosomes dynamically from the BAM file
    autosomal_chromosomes=($(get_autosomal_chromosomes "$bam_file"))

    # Check high-quality reads for each autosomal chromosome
    valid_file=true
    for chr in "${autosomal_chromosomes[@]}"; do
      # Count high-quality reads for each chromosome
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
      # Extract the identifier and add it to the temporary filtered list
      identifier=$(basename "$bam_file" .bam)
      echo "$identifier" >> "$temp_filtered_list"
    fi
  done

  # Count the number of BAM files in the filtered list
  filtered_bam_count=$(wc -l < "$temp_filtered_list")
  echo "Number of BAM files in the filtered list: $filtered_bam_count"

  # If filtering was successful, move the temporary list to the final destination
  if [ "$filtered_bam_count" -gt 0 ]; then
    mv "$temp_filtered_list" "$filtered_list"
    echo "Filtered list of BAM files has been created."
  else
    echo "No valid BAM files found. Not saving filtered list."
    rm "$temp_filtered_list"
  fi

else
  echo "Filtered list already exists. Skipping filtering step."
fi

