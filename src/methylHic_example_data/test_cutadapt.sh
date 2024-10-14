#!/bin/bash

force_download=TRUE 

eval "$(conda shell.bash hook)"
conda activate sra_env

# Define working directory and other important variables
working_dir="/oceanus/collab/InternalJeff/users/kxc732/downloads/projects/methyHic"
SRR="SRR7770845" # Example SRR, replace with appropriate SRR ID if necessary
SRA_FILE="$working_dir/${SRR}.sra"
FASTQ_R1="$working_dir/${SRR}_1.fastq.gz"
FASTQ_R2="$working_dir/${SRR}_2.fastq.gz"
EXTRACTED_R1="$working_dir/${SRR}_extracted_1.fastq.gz"
EXTRACTED_R2="$working_dir/${SRR}_extracted_2.fastq.gz"
SRR_DIR="$working_dir/$SRR"
BARCODE_FILE="$working_dir/barcodes.txt" # Replace with the correct path to your barcode file
DEMUX_DIR="$working_dir/demux" # Directory for demultiplexed files

# Create necessary directories
mkdir -p "$SRR_DIR" "$DEMUX_DIR"

# Remove incomplete downloads or force re-download
if [ "$force_download" = TRUE ]; then
    echo "Force download is enabled. Removing existing SRA and FASTQ files for ${SRR}."
    rm -f "${SRA_FILE}" "${FASTQ_R1}" "${FASTQ_R2}"
    rm -rf "${SRR_DIR}"
fi

# Download SRA file if not present
if [ ! -f "${SRA_FILE}" ]; then
    echo "Downloading ${SRR}..."
    mkdir -p "${SRR_DIR}"
    echo "Running prefetch for ${SRR}..."
    # require 2.10 or greater 
    prefetch --max-size 100000000000 -O "${working_dir}/" "${SRR}" > "${SRR_DIR}/prefetch.log" 2>&1
    PREFETCH_EXIT_CODE=$?
    if [ $PREFETCH_EXIT_CODE -ne 0 ]; then
        echo "Error: Download of ${SRR} failed with exit code $PREFETCH_EXIT_CODE. See ${SRR_DIR}/prefetch.log for details."
        rm -rf "${SRR_DIR}"
        exit 1
    fi
else
    echo "${SRR}.sra already exists in ${SRR_DIR}. Skipping download."
fi

# Convert SRA to FASTQ
if [ ! -f "${FASTQ_R1}" ] || [ ! -f "${FASTQ_R2}" ]; then
    echo "Converting ${SRR}.sra to FASTQ format using fasterq-dump..."
    echo "Running fasterq-dump for ${SRR}..."
    cd "${SRR_DIR}"
    fasterq-dump --split-files "${SRR}" -O "${working_dir}/" > "${SRR_DIR}/fasterq-dump.log" 2>&1
    FASTERQ_EXIT_CODE=$?
    cd -  # Return to the previous directory
    if [ $FASTERQ_EXIT_CODE -ne 0 ]; then
        echo "Error: FASTQ conversion for ${SRR}.sra failed with exit code $FASTERQ_EXIT_CODE. See ${SRR_DIR}/fasterq-dump.log for details."
        rm -f "${working_dir}/${SRR}_1.fastq" "${working_dir}/${SRR}_2.fastq"
        exit 1
    fi

    # Verify FASTQ files are created and non-empty
    if [ ! -s "${working_dir}/${SRR}_1.fastq" ] || [ ! -s "${working_dir}/${SRR}_2.fastq" ]; then
        echo "Error: FASTQ files for ${SRR} are empty or not created."
        rm -f "${working_dir}/${SRR}_1.fastq" "${working_dir}/${SRR}_2.fastq"
        exit 1
    fi

    echo "Compressing FASTQ files for ${SRR}..."
    gzip "${working_dir}/${SRR}_1.fastq" "${working_dir}/${SRR}_2.fastq"
else
    echo "FASTQ files for ${SRR} already exist. Skipping conversion."
fi

# # Verify FASTQ files are valid
# echo "Verifying integrity of FASTQ files for ${SRR}..."
# if ! gunzip -t "${FASTQ_R1}" || ! gunzip -t "${FASTQ_R2}"; then
#     echo "Error: FASTQ files for ${SRR} are corrupted. Deleting and retrying."
#     rm -f "${FASTQ_R1}" "${FASTQ_R2}"
#     exit 1
# fi

# # Extract UMI and cell barcodes
# if [ ! -f "${EXTRACTED_R1}" ] || [ ! -f "${EXTRACTED_R2}" ]; then
#     echo "Extracting UMI and cell barcodes for ${SRR}..."
#     echo "Running umi_tools extract for ${SRR}..."
#     umi_tools extract \
#         --stdin "${FASTQ_R1}" \
#         --read2-in="${FASTQ_R2}" \
#         --bc-pattern=NNNNNNNNNNCCCCCCCCC \
#         --stdout="${EXTRACTED_R1}" \
#         --read2-out="${EXTRACTED_R2}" > "${SRR_DIR}/umi_tools_extract.log" 2>&1
#     UMI_EXTRACT_EXIT_CODE=$?
#     if [ $UMI_EXTRACT_EXIT_CODE -ne 0 ]; then
#         echo "Error: UMI extraction for ${SRR} failed with exit code $UMI_EXTRACT_EXIT_CODE. See ${SRR_DIR}/umi_tools_extract.log for details."
#         rm -f "${EXTRACTED_R1}" "${EXTRACTED_R2}"
#         exit 1
#     fi

#     # Verify extracted files are created and non-empty
#     if [ ! -s "${EXTRACTED_R1}" ] || [ ! -s "${EXTRACTED_R2}" ]; then
#         echo "Error: Extracted FASTQ files for ${SRR} are empty or not created."
#         rm -f "${EXTRACTED_R1}" "${EXTRACTED_R2}"
#         exit 1
#     fi
# else
#     echo "Extracted UMI and cell barcodes for ${SRR} already exist. Skipping extraction."
# fi

# Demultiplex using cutadapt
echo "Demultiplexing FASTQ files for ${SRR} using cutadapt..."
echo "Running cutadapt for ${SRR}..."
cutadapt -g "file:${BARCODE_FILE}" \
    -o "${DEMUX_DIR}/${SRR}_demux_{name}_R1.fastq.gz" \
    -p "${DEMUX_DIR}/${SRR}_demux_{name}_R2.fastq.gz" \
    "${EXTRACTED_R1}" "${EXTRACTED_R2}" > "${SRR_DIR}/cutadapt.log" 2>&1
CUTADAPT_EXIT_CODE=$?
if [ $CUTADAPT_EXIT_CODE -ne 0 ]; then
    echo "Error: Demultiplexing for ${SRR} failed with exit code $CUTADAPT_EXIT_CODE. See ${SRR_DIR}/cutadapt.log for details."
    exit 1
fi

# echo "Demultiplexing for ${SRR} completed successfully."
