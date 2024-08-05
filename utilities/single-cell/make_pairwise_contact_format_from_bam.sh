#!/bin/bash

eval "$(python3 config_and_print.py)"

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

# Process each BAM file in the output directory
for bam_file in $output_directory/sc*.bam; do
  echo "Processing $bam_file"
  
  # Extract prefix from filename for job naming and command construction
  prefix=$(basename "$bam_file" .bam)

  # Construct the command to be executed
  good_reads_bam="$output_directory/${prefix}.good_reads.bam"
  hic_txt="$output_directory/${prefix}.hic.txt"

  # Skip processing if output files already exist
  if [ -f "$good_reads_bam" ] && [ -f "$hic_txt" ]; then
    echo "Output files for $bam_file already exist. Skipping."
    continue
  fi

  # Filter BAM file and convert to juicer_short format
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
done

echo "All BAM files have been processed."

