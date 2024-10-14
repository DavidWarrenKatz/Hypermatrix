#!/bin/bash

# Load necessary modules
module load sra-toolkit/2.10.9

# Create necessary directories
mkdir -p ../../projects/methyHic
mkdir -p ../../projects/methyHic/demux  # Directory for demultiplexed FASTQ files
mkdir -p ../../projects/methyHic/extracted  # Directory for extracted FASTQ files
 
# Create a log file for failed cutadapt processes
FAILED_LOG="cutadapt_failed.log"
> $FAILED_LOG  # Clear the file if it exists

# Function to download, convert, extract UMI/barcodes, and demultiplex SRA to FASTQ
process_srr() {
    SRR=$1
    echo "[$(date)] Processing ${SRR} with PID $$"

    # Check if the directory for the SRR ID exists and skip download if it does
    if [ ! -d "../../projects/methyHic/${SRR}" ]; then
        echo "Downloading ${SRR}..."
        prefetch --max-size 100000000000 -O ../../projects/methyHic/ ${SRR}
        if [ $? -ne 0 ]; then
            echo "Error: Download of ${SRR} failed."
            return 1
        fi
    else
        echo "${SRR}.sra already exists in ${SRR} directory. Skipping download."
    fi

    # Convert SRA to FASTQ if not already converted
    if [ ! -f "../../projects/methyHic/${SRR}_1.fastq.gz" ] && [ ! -f "../../projects/methyHic/${SRR}_2.fastq.gz" ]; then
        echo "Converting ${SRR}.sra to FASTQ format..."
        fastq-dump --split-files --gzip "../../projects/methyHic/${SRR}/${SRR}.sra" -O ../../projects/methyHic/
        if [ $? -ne 0 ]; then
            echo "Error: FASTQ conversion for ${SRR}.sra failed."
            return 1
        fi
    else
        echo "FASTQ files for ${SRR} already exist. Skipping conversion."
    fi

    # Extract UMI and cell barcodes from Read 1
    if [ ! -f "../../projects/methyHic/extracted/${SRR}_R1_extracted.fastq.gz" ] || [ ! -f "../../projects/methyHic/extracted/${SRR}_R2_extracted.fastq.gz" ]; then
        echo "Extracting UMI and cell barcodes for ${SRR}..."
        umi_tools extract \
            --stdin ../../projects/methyHic/${SRR}_1.fastq.gz \
            --read2-in=../../projects/methyHic/${SRR}_2.fastq.gz \
            --bc-pattern=NNNNNNNNNNCCCCCCCCC \
            --stdout=../../projects/methyHic/extracted/${SRR}_R1_extracted.fastq.gz \
            --read2-out=../../projects/methyHic/extracted/${SRR}_R2_extracted.fastq.gz

        if [ $? -ne 0 ]; then
            echo "Error: UMI extraction for ${SRR} failed."
            return 1
        fi
    else
        echo "Extracted UMI and cell barcodes for ${SRR} already exist. Skipping extraction."
    fi

    # Define barcode sequences for demultiplexing (AD002, AD006, AD008, AD010)
    BARCODE_FILE=barcodes.fasta

    # Check if demultiplexed files already exist before running cutadapt
    if [ ! -f "../../projects/methyHic/demux/${SRR}_demux_AD002_R1.fastq.gz" ] || [ ! -f "../../projects/methyHic/demux/${SRR}_demux_AD002_R2.fastq.gz" ]; then
   
        # Run the cutadapt command and handle errors
        cutadapt -g file:${BARCODE_FILE} \
            -o ../../projects/methyHic/demux/${SRR}_demux_{name}_R1.fastq.gz \
            -p ../../projects/methyHic/demux/${SRR}_demux_{name}_R2.fastq.gz \
            ../../projects/methyHic/extracted/${SRR}_R1_extracted.fastq.gz \
            ../../projects/methyHic/extracted/${SRR}_R2_extracted.fastq.gz

        if [ $? -ne 0 ]; then
            echo "Error: Demultiplexing for ${SRR} failed."
            echo "${SRR}" >> $FAILED_LOG  # Log the failed SRR ID
            return 1
        fi
    else
        echo "Demultiplexed FASTQ files for ${SRR} already exist. Skipping cutadapt."
    fi

    echo "Demultiplexing for ${SRR} completed."
}

export -f process_srr

# Use GNU Parallel to download, convert, extract UMI/barcodes, and demultiplex SRA files in parallel
cat SRR_Acc_List_sc.txt | parallel -j 8 process_srr

echo "All downloads, conversions, UMI extractions, and demultiplexing completed."

# Check for failed SRR IDs
if [ -s $FAILED_LOG ]; then
    echo "The following SRR IDs failed during demultiplexing:"
    cat $FAILED_LOG
else
    echo "All SRR IDs processed successfully."
fi

