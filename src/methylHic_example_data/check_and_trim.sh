#!/bin/bash

module load FastQC/0.11.9
module load Trimmomatic/0.39

# Check for required tools
if ! command -v fastqc &> /dev/null || ! command -v trimmomatic &> /dev/null
then
    echo "FastQC and/or Trimmomatic not found. Please install them and try again."
    exit 1
fi

# Directory containing FASTQ files
fastq_dir="/jet/home/dkatz/tmp_ondemand_ocean_mcb190124p_symlink/dkatz/hypermatrix/projects/methyHic"

# Change to the directory containing FASTQ files
cd "$fastq_dir" || { echo "Directory not found: $fastq_dir"; exit 1; }

# Input argument: the name of the example FASTQ file
example_fastq="$fastq_dir/SRR7770812_1.fastq.gz"

# Check if the example FASTQ file is provided
if [ ! -f "$example_fastq" ]; then
  echo "Example FASTQ file not found: $example_fastq"
  exit 1
fi

# Run FastQC on the example FASTQ file
echo "Running FastQC on example file: $example_fastq"
fastqc -o . --nogroup "$example_fastq"

# Check the FastQC report for adapter contamination
report_file="${example_fastq%.fastq.gz}_fastqc/fastqc_data.txt"
if grep -q "Adapter Content" "$report_file"; then
  adapter_content=$(grep "Adapter Content" "$report_file" | tail -n 1 | awk '{print $NF}')
  echo "Adapter content detected: $adapter_content%"
  if (( $(echo "$adapter_content > 5" | bc -l) )); then
    echo "Adapters detected. The reads have not been trimmed. Proceeding with trimming."

    # Loop over all paired-end FASTQ files and trim them using Trimmomatic
    for fastq_file_1 in *_1.fastq.gz; do
      fastq_file_2="${fastq_file_1/_1.fastq.gz/_2.fastq.gz}"
      
      # Check if the paired-end file exists
      if [ ! -f "$fastq_file_2" ]; then
        echo "Paired file not found for $fastq_file_1"
        continue
      fi

      echo "Processing $fastq_file_1 and $fastq_file_2"
      
      # Define output files
      output_prefix="${fastq_file_1%_1.fastq.gz}"
      paired_output1="${output_prefix}_paired_R1.fastq.gz"
      unpaired_output1="${output_prefix}_unpaired_R1.fastq.gz"
      paired_output2="${output_prefix}_paired_R2.fastq.gz"
      unpaired_output2="${output_prefix}_unpaired_R2.fastq.gz"
      
      # Run Trimmomatic with default parameters for paired-end trimming
      trimmomatic PE \
        -threads 4 \
        "$fastq_file_1" "$fastq_file_2" \
        "$paired_output1" "$unpaired_output1" \
        "$paired_output2" "$unpaired_output2" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
      
      echo "Trimming completed for $fastq_file_1 and $fastq_file_2"
    done
    
  else
    echo "No significant adapter contamination detected. The reads appear to be already trimmed."
  fi
else
  echo "Error: Could not determine adapter content from FastQC report."
  exit 1
fi

