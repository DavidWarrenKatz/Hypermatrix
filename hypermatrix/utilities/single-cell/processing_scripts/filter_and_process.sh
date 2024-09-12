######################################################################################
# This shell script iterates through the bam_directory and makes a filtered
# list of bam files that have a sufficient number of high quality reads.
# It also filters the bam files to only contain good reads and processes them into juicer format
# with the function sam2juicer.py. Symbolic links of all the bam files that pass filtering are made into the output directory. 
# If the filtered list is already present in the output directory, this step is skipped.  
# Author: David Katz (davidkatz02@gmail.com)                                                              
#####################################################################################

#!/bin/bash

eval "$(python3 config_and_print.py)"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Source the conda environment setup script
source /software/miniconda3/4.12.0/etc/profile.d/conda.sh

# Define environment names
schicluster_env=schicluster2
bisulfite_env=bisulfitehic27

# Activate the desired conda environment
conda activate $bisulfite_env

# Load samtools module
module load samtools

# Ensure that bisulfite_env_python_path is set correctly
bisulfite_env_python_path=$(which python)

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

      # Create symbolic link
      if [ ! -L "$output_directory/$identifier.bam" ]; then
        ln -s "$bam_directory/$identifier.b37.calmd.bam" "$output_directory/$identifier.bam"
      fi

      # Process BAM file and convert to juicer_short format
      good_reads_bam="$output_directory/${identifier}.good_reads.bam"
      hic_txt="$output_directory/${identifier}.hic.txt"

      # Skip processing if output files already exist
      if [ -f "$good_reads_bam" ] && [ -f "$hic_txt" ]; then
        echo "Output files for $bam_file already exist. Skipping."
        continue
      fi

      # Filter BAM file
      samtools view --threads 5 -bh -q 30 -f 1 -F 1804 "$bam_file" > "$good_reads_bam"
      if [ $? -ne 0 ]; then
        echo "ERROR processing $bam_file with samtools."
        continue
      fi

      # Run sam2juicer.py script
      "$bisulfite_env_python_path" "$software_directory/sam2juicer.py" -s "$good_reads_bam" -f "$fragments_file" > "$hic_txt"
      if [ $? -ne 0 ]; then
        echo "WARNING: Problem processing $good_reads_bam with sam2juicer.py. Check the fragments file: $fragments_file"
        continue
      fi

      echo "Processed $bam_file -> $hic_txt"
    fi
  done

  # Count the number of BAM files in the filtered list
  filtered_bam_count=$(wc -l < "$filtered_list")
  echo "Number of BAM files in the filtered list: $filtered_bam_count"
  echo "Filtered list of BAM files has been created."
else
  echo "Filtered list already exists. Skipping filtering step."
fi

echo "All BAM files have been processed."

