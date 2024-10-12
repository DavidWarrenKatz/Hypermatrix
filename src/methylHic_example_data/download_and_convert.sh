#!/bin/bash
# file: download_and_convert.sh

# Enable debugging to trace variables and commands
set -x

# Force download variable
force_download=TRUE

# Function to check if a package is installed in the environment
function check_and_install_package() {
    package_name=$1
    required_version=$2
    channel=$3

    echo "Checking if $package_name version $required_version is installed..."
    if ! conda list | grep -q "^$package_name\s*$required_version"; then
        echo "$package_name version $required_version not found. Installing..."
        conda install -y -c $channel $package_name=$required_version || { echo "Failed to install $package_name $required_version"; exit 1; }
    else
        echo "$package_name version $required_version is already installed."
    fi
}

# Check if conda is available and install necessary packages
echo "Checking if conda is available..."
if ! command -v conda &> /dev/null; then
    echo "Conda not found! Please install Miniconda or Anaconda first."
    exit 1
fi

# Create a new conda environment if it doesn't exist
echo "Checking if conda environment 'sra_env' exists..."
if ! conda env list | grep -q 'sra_env'; then
    echo "Creating conda environment 'sra_env' and installing necessary tools..."
    conda create -n sra_env -y python=3.8 || { echo "Failed to create conda environment and install tools"; exit 1; }
else
    echo "Conda environment 'sra_env' already exists."
fi

# Activate the 'sra_env' environment and install necessary packages
echo "Activating conda environment 'sra_env'..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra_env

# Check and install necessary packages if not already installed
echo "Ensuring necessary packages are installed..."
check_and_install_package "sra-tools" "2.10.7" "bioconda"
check_and_install_package "umi_tools" "1.1.5" "bioconda"
check_and_install_package "cutadapt" "4.9" "bioconda"

# Base directories
echo "Setting up base directories..."
BASE_DIR="/oceanus/collab/InternalJeff/users/kxc732/downloads/projects/methyHic"
DEMUX_DIR="${BASE_DIR}/demux"
EXTRACTED_DIR="${BASE_DIR}/extracted"

echo "Creating necessary directories..."
mkdir -p "${DEMUX_DIR}"
mkdir -p "${EXTRACTED_DIR}"

# Define and export barcode file path
BARCODE_FILE="/oceanus/collab/InternalJeff/users/kxc732/downloads/new/Hypermatrix/src/methylHic_example_data/barcodes.fasta"
export BARCODE_FILE

echo "Checking if barcode file exists at $BARCODE_FILE..."
if [ ! -f "$BARCODE_FILE" ]; then
    echo "Error: Barcode file $BARCODE_FILE not found!"
    exit 1
fi
echo "Using barcode file: $BARCODE_FILE"

# Export variables needed for parallel processing
export BASE_DIR DEMUX_DIR EXTRACTED_DIR BARCODE_FILE force_download

# Function to process SRR IDs
process_srr() {
    SRR=$1
    echo "[$(date)] Processing ${SRR} with PID $$"

    SRR_DIR="${BASE_DIR}/${SRR}"
    SRA_FILE="${SRR_DIR}/${SRR}.sra"
    FASTQ_R1="${BASE_DIR}/${SRR}_1.fastq.gz"
    FASTQ_R2="${BASE_DIR}/${SRR}_2.fastq.gz"
    EXTRACTED_R1="${EXTRACTED_DIR}/${SRR}_R1_extracted.fastq.gz"
    EXTRACTED_R2="${EXTRACTED_DIR}/${SRR}_R2_extracted.fastq.gz"

    echo "SRR_DIR: $SRR_DIR"
    echo "SRA_FILE: $SRA_FILE"
    echo "FASTQ_R1: $FASTQ_R1"
    echo "FASTQ_R2: $FASTQ_R2"
    echo "EXTRACTED_R1: $EXTRACTED_R1"
    echo "EXTRACTED_R2: $EXTRACTED_R2"

    # Remove incomplete downloads or force re-download
    if [ "$force_download" = TRUE ]; then
        echo "Force download is enabled. Removing existing SRA and FASTQ files for ${SRR}."
        rm -f "${SRA_FILE}"
        rm -f "${FASTQ_R1}" "${FASTQ_R2}"
        rm -rf "${SRR_DIR}"
    elif [ -d "${SRR_DIR}" ] && [ ! -f "${SRA_FILE}" ]; then
        echo "Incomplete download detected for ${SRR}. Removing directory and retrying download."
        rm -rf "${SRR_DIR}"
    fi

    # Download SRA file if not present
    if [ ! -f "${SRA_FILE}" ]; then
        echo "Downloading ${SRR}..."
        mkdir -p "${SRR_DIR}"
        echo "Running prefetch for ${SRR}..."
        prefetch --max-size 100000000000 -O "${BASE_DIR}/" "${SRR}" > "${SRR_DIR}/prefetch.log" 2>&1
        PREFETCH_EXIT_CODE=$?
        if [ $PREFETCH_EXIT_CODE -ne 0 ]; then
            echo "Error: Download of ${SRR} failed with exit code $PREFETCH_EXIT_CODE. See ${SRR_DIR}/prefetch.log for details."
            rm -rf "${SRR_DIR}"
            return 1
        fi
    else
        echo "${SRR}.sra already exists in ${SRR_DIR}. Skipping download."
    fi

    # Verify SRA file exists
    if [ ! -f "${SRA_FILE}" ]; then
        echo "Error: SRA file ${SRA_FILE} does not exist after download."
        return 1
    fi

    # Convert SRA to FASTQ
    if [ ! -f "${FASTQ_R1}" ] || [ ! -f "${FASTQ_R2}" ]; then
        echo "Converting ${SRR}.sra to FASTQ format using fasterq-dump..."
        echo "Running fasterq-dump for ${SRR}..."
        cd "${SRR_DIR}"
        fasterq-dump --split-files "${SRR}" -O "${BASE_DIR}/" > "${SRR_DIR}/fasterq-dump.log" 2>&1
        FASTERQ_EXIT_CODE=$?
        cd -  # Return to the previous directory
        if [ $FASTERQ_EXIT_CODE -ne 0 ]; then
            echo "Error: FASTQ conversion for ${SRR}.sra failed with exit code $FASTERQ_EXIT_CODE. See ${SRR_DIR}/fasterq-dump.log for details."
            rm -f "${BASE_DIR}/${SRR}_1.fastq" "${BASE_DIR}/${SRR}_2.fastq"
            return 1
        fi

        # Verify FASTQ files are created and non-empty
        if [ ! -s "${BASE_DIR}/${SRR}_1.fastq" ] || [ ! -s "${BASE_DIR}/${SRR}_2.fastq" ]; then
            echo "Error: FASTQ files for ${SRR} are empty or not created."
            rm -f "${BASE_DIR}/${SRR}_1.fastq" "${BASE_DIR}/${SRR}_2.fastq"
            return 1
        fi

        echo "Compressing FASTQ files for ${SRR}..."
        gzip "${BASE_DIR}/${SRR}_1.fastq" "${BASE_DIR}/${SRR}_2.fastq"
    else
        echo "FASTQ files for ${SRR} already exist. Skipping conversion."
    fi

    # Verify FASTQ files are valid
    echo "Verifying integrity of FASTQ files for ${SRR}..."
    if ! gunzip -t "${FASTQ_R1}" || ! gunzip -t "${FASTQ_R2}"; then
        echo "Error: FASTQ files for ${SRR} are corrupted. Deleting and retrying."
        rm -f "${FASTQ_R1}" "${FASTQ_R2}"
        return 1
    fi

    # Extract UMI and cell barcodes
    if [ ! -f "${EXTRACTED_R1}" ] || [ ! -f "${EXTRACTED_R2}" ]; then
        echo "Extracting UMI and cell barcodes for ${SRR}..."
        echo "Running umi_tools extract for ${SRR}..."
        umi_tools extract \
            --stdin "${FASTQ_R1}" \
            --read2-in="${FASTQ_R2}" \
            --bc-pattern=NNNNNNNNNNCCCCCCCCC \
            --stdout="${EXTRACTED_R1}" \
            --read2-out="${EXTRACTED_R2}" > "${SRR_DIR}/umi_tools_extract.log" 2>&1
        UMI_EXTRACT_EXIT_CODE=$?
        if [ $UMI_EXTRACT_EXIT_CODE -ne 0 ]; then
            echo "Error: UMI extraction for ${SRR} failed with exit code $UMI_EXTRACT_EXIT_CODE. See ${SRR_DIR}/umi_tools_extract.log for details."
            rm -f "${EXTRACTED_R1}" "${EXTRACTED_R2}"
            return 1
        fi

        # Verify extracted files are created and non-empty
        if [ ! -s "${EXTRACTED_R1}" ] || [ ! -s "${EXTRACTED_R2}" ]; then
            echo "Error: Extracted FASTQ files for ${SRR} are empty or not created."
            rm -f "${EXTRACTED_R1}" "${EXTRACTED_R2}"
            return 1
        fi
    else
        echo "Extracted UMI and cell barcodes for ${SRR} already exist. Skipping extraction."
    fi

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
        return 1
    fi

    echo "Demultiplexing for ${SRR} completed successfully."
}

export -f process_srr

# Process SRR IDs using GNU Parallel
echo "Starting parallel processing of SRR IDs..."
cat SRR_Acc_List_sc.txt | parallel -j 64 process_srr

echo "All downloads, conversions, UMI extractions, and demultiplexing completed."

# Disable debugging
set +x
