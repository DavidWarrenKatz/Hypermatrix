#./make_pairwise_contact_format_from_bam.sh

######################################################################################                   #This shell script iterates through all the symbolic links in the output directory.
#The bam files are filtered to only contain good reads, and then processed into juicer format
#with the function sam2juicer.py
#Author: David Katz (davidkatz02@gmail.com)                                                              #####################################################################################

#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Import the parameters from config.py (relative to the script's directory)
eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

# Load samtools module
module load samtools

# Process each BAM file in the output directory
for bam_file in $output_directory/sc*.bam; do
  # Skip files containing ".good_reads.bam"
  if [[ "$bam_file" == *".good_reads.bam" ]]; then
    continue
  fi
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
  python "$software_directory/sam2juicer.py" -s "$good_reads_bam" -f "$fragments_file" > "$hic_txt"
  if [ $? -ne 0 ]; then
    echo "WARNING: Problem processing $good_reads_bam with sam2juicer.py. Check the fragments file: $fragments_file"
    continue
  fi
  
  echo "Processed $bam_file -> $hic_txt"
done

echo "All BAM files have been processed."

