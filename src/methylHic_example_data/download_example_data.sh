#!/bin/bash

#Do this if parallel not available
# Download and install GNU Parallel locally
#wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
#tar -xjf parallel-latest.tar.bz2
#cd parallel-*
#./configure --prefix=$HOME/.local
#make
#make install

# Load necessary modules
module load sra-toolkit/2.10.9
module load FastQC/0.11.9

# Create directory if it doesn't exist
mkdir -p ../../projects/methyHic
mkdir -p ../../projects/methyHic/fastqc_reports

# Check if the main tar file is already downloaded
if [ ! -f ../../projects/methyHic/GSE119171_RAW.tar ]; then
    echo "Downloading GSE119171_RAW.tar..."
    wget -O ../../projects/methyHic/GSE119171_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119171/suppl/GSE119171_RAW.tar"
else
    echo "GSE119171_RAW.tar already exists. Skipping download."
fi

# Check if the files from the tar archive are already extracted
if [ ! -f "../../projects/methyHic/GSM3359999_JL_457_1.CGATGT.mm9.calmd.cpg.filtered.sort.CG.strand.6plus2.bed.gz" ]; then
    echo "Extracting GSE119171_RAW.tar..."
    tar -xvf ../../projects/methyHic/GSE119171_RAW.tar -C ../../projects/methyHic/
else
    echo "Files already extracted. Skipping extraction."
fi

# Function to download, convert, and verify each SRR
process_srr() {
    SRR=$1
    echo "[$(date)] Processing ${SRR} with PID $$"
    
    # Check if the directory for the SRR ID exists and skip download if it does
    if [ ! -d "../../projects/methyHic/${SRR}" ]; then
        echo "Downloading ${SRR}..."
        prefetch -O ../../projects/methyHic/ ${SRR}
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

    # Run FastQC to check the integrity of the FASTQ files
    echo "Running FastQC on ${SRR} FASTQ files..."
    fastqc ../../projects/methyHic/${SRR}_*.fastq.gz -o ../../projects/methyHic/fastqc_reports/
    if [ $? -ne 0 ]; then
        echo "Error: FastQC failed for ${SRR}."
        return 1
    fi

    echo "[$(date)] Finished processing ${SRR} with PID $$"
}

export -f process_srr

# Use GNU Parallel to download, convert, and verify SRA files in parallel
# Logging the start and end of each job
cat SRR_Acc_List.txt | parallel -j 8 process_srr

echo "All downloads, conversions, and verifications completed."

