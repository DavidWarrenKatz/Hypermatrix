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

# Barcode file path
BARCODE_FILE="../../projects/methyHic/barcodes.fasta"

# Create the barcode file in FASTA format
echo "Creating barcode file in FASTA format..."
cat <<EOL > $BARCODE_FILE
>sample1
CGATGT
>sample2
GCCAAT
>sample3
ACTTGA
>sample4
TAGCTT
EOL

echo "Barcode file created at $BARCODE_FILE"

# Function to download, convert, and demultiplex SRA to FASTQ
process_srr() {
    SRR=$1
    echo "[$(date)] Processing ${SRR} with PID $$"

    SRR_DIR="../../projects/methyHic/${SRR}"
    SRA_FILE="${SRR_DIR}/${SRR}.sra"
    FASTQ_R1="../../projects/methyHic/${SRR}_1.fastq.gz"
    FASTQ_R2="../../projects/methyHic/${SRR}_2.fastq.gz"
    DEMUX_DIR="../../projects/methyHic/demux"

    # Barcode file path
    BARCODE_FILE="../../projects/methyHic/barcodes.fasta"

    # Check if the directory for the SRR ID exists and skip download if it does
    if [ ! -d "$SRR_DIR" ]; then
        echo "Downloading ${SRR}..."
        prefetch --max-size 100000000000 -O ../../projects/methyHic/ ${SRR}
        if [ $? -ne 0 ]; then
            echo "Error: Download of ${SRR} failed."
            return 1
        fi
    else
        echo "${SRR}.sra already exists in ${SRR} directory. Skipping download."
    fi

    # Convert SRA to FASTQ
    if [ ! -f "${FASTQ_R1}" ] || [ ! -f "${FASTQ_R2}" ]; then
        echo "Converting ${SRR}.sra to FASTQ format using fasterq-dump..."
        fasterq-dump --split-files --gzip -O ../../projects/methyHic/ ${SRR}
        if [ $? -ne 0 ]; then
            echo "Error: Conversion of ${SRR}.sra to FASTQ failed."
            return 1
        fi
    else
        echo "FASTQ files for ${SRR} already exist. Skipping conversion."
    fi

    # Demultiplex using cutadapt
    echo "Demultiplexing FASTQ files for ${SRR} using cutadapt..."
    cutadapt -g "file:${BARCODE_FILE}" \
        -o "${DEMUX_DIR}/${SRR}_demux_{name}_R1.fastq.gz" \
        -p "${DEMUX_DIR}/${SRR}_demux_{name}_R2.fastq.gz" \
        "${FASTQ_R1}" "${FASTQ_R2}" > "${SRR_DIR}/cutadapt.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: Demultiplexing for ${SRR} failed. Check ${SRR_DIR}/cutadapt.log for details."
        echo "${SRR}" >> $FAILED_LOG
        return 1
    fi

    echo "Demultiplexing for ${SRR} completed successfully."
}

export -f process_srr

# Process SRR IDs from the list using GNU Parallel
echo "Starting parallel processing of SRR IDs..."
cat SRR_Acc_List_sc.txt | parallel -j 64 process_srr

echo "All downloads, conversions, and demultiplexing completed."

