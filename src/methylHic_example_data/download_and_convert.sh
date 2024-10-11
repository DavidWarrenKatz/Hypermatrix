#!/bin/bash

# Load necessary modules
module load sra-toolkit/2.10.9

# Create necessary directories
mkdir -p ../../projects/methyHic
mkdir -p ../../projects/methyHic/demux  # Directory for demultiplexed FASTQ files


# Define barcode sequences for demultiplexing (AD002, AD006, AD008, AD010)
BARCODE_FILE=barcodes.fasta

# Create the barcode file
cat > ${BARCODE_FILE} << EOL
>P5L_AD002
TTCCCTACACGACGCTCTTCCGATCTCGATGT
>P5L_AD006
TTCCCTACACGACGCTCTTCCGATCTGCCAAT
>P5L_AD008
TTCCCTACACGACGCTCTTCCGATCTACTTGA
>P5L_AD010
TTCCCTACACGACGCTCTTCCGATCTTAGCTT
EOL

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
    echo "Extracting UMI and cell barcodes for ${SRR}..."
    umi_tools extract \
        --stdin ../../projects/methyHic/${SRR}_1.fastq.gz \
        --read2-in=../../projects/methyHic/${SRR}_2.fastq.gz \
        --bc-pattern=NNNNNNNNNNCCCCCCCCC \
        --stdout=../../projects/methyHic/demux/${SRR}_R1_extracted.fastq.gz \
        --read2-out=../../projects/methyHic/demux/${SRR}_R2_extracted.fastq.gz

    if [ $? -ne 0 ]; then
        echo "Error: UMI extraction for ${SRR} failed."
        return 1
    fi

    # Demultiplex the extracted FASTQ files based on cell barcodes
    echo "Demultiplexing FASTQ files for ${SRR}..."
    cutadapt -g file:${BARCODE_FILE} \
        -o ../../projects/methyHic/demux/${SRR}_demux_{name}_R1.fastq.gz \
        -p ../../projects/methyHic/demux/${SRR}_demux_{name}_R2.fastq.gz \
        ../../projects/methyHic/demux/${SRR}_R1_extracted.fastq.gz \
        ../../projects/methyHic/demux/${SRR}_R2_extracted.fastq.gz

    if [ $? -ne 0 ]; then
        echo "Error: Demultiplexing for ${SRR} failed."
        return 1
    fi

    echo "Demultiplexing for ${SRR} completed."
}

export -f process_srr

# Use GNU Parallel to download, convert, extract UMI/barcodes, and demultiplex SRA files in parallel
cat SRR_Acc_List.txt | parallel -j 8 process_srr

echo "All downloads, conversions, UMI extractions, and demultiplexing completed."

