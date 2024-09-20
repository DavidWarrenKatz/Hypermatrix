#!/bin/bash

module load FastQC/0.11.9
module load Trimmomatic/0.39

# Exit on any error
set -e

# Directory containing FASTQ files
fastq_dir="/jet/home/dkatz/tmp_ondemand_ocean_mcb190124p_symlink/dkatz/hypermatrix/projects/methyHic"

# Change to the directory containing FASTQ files
cd "$fastq_dir" || { echo "Directory not found: $fastq_dir"; exit 1; }

# Loop over all paired-end FASTQ files
for fastq_file_1 in *_1.fastq.gz; do
    fastq_file_2="${fastq_file_1/_1.fastq.gz/_2.fastq.gz}"
    
    # Check if the paired-end file exists
    if [ ! -f "$fastq_file_2" ]; then
        echo "Paired file not found for $fastq_file_1"
        continue
    fi

    output_prefix="${fastq_file_1%_1.fastq.gz}"
    fastqc_zip="${fastq_dir}/${output_prefix}_fastqc.zip"
    fastqc_report="${fastq_dir}/${output_prefix}_fastqc/fastqc_data.txt"
    paired_output1="${output_prefix}_paired_R1.fastq.gz"
    paired_output2="${output_prefix}_paired_R2.fastq.gz"
    
    # Check if FastQC was already run by verifying both the zip and the report exist
    if [ -f "$fastqc_report" ]; then
        echo "FastQC report already exists for $fastq_file_1. Skipping FastQC..."
    elif [ -f "$fastqc_zip" ]; then
        echo "Unzipping existing FastQC zip file: $fastqc_zip"
        unzip -o "$fastqc_zip"
    else
        echo "Running FastQC on $fastq_file_1"
        fastqc -o "$fastq_dir" --nogroup "$fastq_file_1"
        
        # Verify FastQC ran successfully by checking the zip file
        if [ ! -f "$fastqc_zip" ]; then
            echo "FastQC zip file not found for $fastq_file_1"
            exit 1
        fi
        
        echo "Unzipping FastQC report: $fastqc_zip"
        unzip -o "$fastqc_zip"
    fi

    # Check if the FastQC report exists after unzipping
    if [ ! -f "$fastqc_report" ]; then
        echo "FastQC report not found for $fastq_file_1. Skipping further analysis."
        continue
    fi

    # Look for adapter content in the FastQC report
    adapter_content=$(grep "Adapter Content" "$fastqc_report" | tail -n 1 | awk '{print $NF}')
    
    # Check if adapter trimming is required and if the trimmed files already exist
    if [ -n "$adapter_content" ] && (( $(echo "$adapter_content > 5" | bc -l) )); then
        echo "Adapter contamination detected: $adapter_content%"

        if [ -f "$paired_output1" ] && [ -f "$paired_output2" ]; then
            echo "Trimmed files already exist for $fastq_file_1. Skipping trimming..."
        else
            echo "Proceeding with trimming $fastq_file_1 and $fastq_file_2"

            unpaired_output1="${output_prefix}_unpaired_R1.fastq.gz"
            unpaired_output2="${output_prefix}_unpaired_R2.fastq.gz"
            
            # Run Trimmomatic for paired-end trimming
            trimmomatic PE \
                -threads 4 \
                "$fastq_file_1" "$fastq_file_2" \
                "$paired_output1" "$unpaired_output1" \
                "$paired_output2" "$unpaired_output2" \
                ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

            echo "Trimming completed for $fastq_file_1 and $fastq_file_2"
        fi
    else
        echo "No significant adapter contamination detected for $fastq_file_1."
    fi
done

