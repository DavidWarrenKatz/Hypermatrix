#!/bin/bash

# Load necessary modules
module load sra-toolkit/2.10.9

# Create necessary directories
mkdir -p ../../projects/methyHic
mkdir -p ../../projects/methyHic/output_bam  # Directory for output BAM files

## Function to download, convert, verify each SRR, and run bisulfite Hi-C mapping
process_srr() {
    SRR=$1
    echo "[$(date)] Processing ${SRR} with PID $$"
    
    # Check if the directory for the SRR ID exists and skip download if it does
    if [ ! -d "../../projects/methyHic/${SRR}" ]; then
        echo "Downloading ${SRR}..."
        prefetch --max-size 100000000000 -O ../../projects/methyHic/ ${SRR}
        if [ $? -ne 0 ]; then
            echo "Error: Download of ${SRR} failed."
            echo "${SRR}" >> $FAILED_LIST
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
            echo "${SRR}" >> $FAILED_LIST
            return 1
        fi
    else
        echo "FASTQ files for ${SRR} already exist. Skipping conversion."
    fi

    # Check if output BAM file already exists
    OUTPUT_BAM="../../projects/methyHic/output_bam/${SRR}_output.bam"
    if [ ! -f "${OUTPUT_BAM}" ]; then
        # Run the bisulfite Hi-C mapping using Java
        echo "Running bisulfite Hi-C mapping for ${SRR}..."

        # Set variables for the paths to make the command more readable
        BISULFITEHIC_DIR="./bisulfitehic/dnaase-bisulfitehic-04680506dd40"
        REFERENCE_GENOME="../../projects/methyHic/mm9.fa"
        FASTQ_R1="../../projects/methyHic/${SRR}_1.fastq.gz"
        FASTQ_R2="../../projects/methyHic/${SRR}_2.fastq.gz"
        ENZYME_LIST="../../projects/methyHic/mm9_DpnII.bed"

        # Java command for bisulfite Hi-C mapping
        java -Xmx15G \
             -Djava.library.path=${BISULFITEHIC_DIR}/jbwa/jbwa-1.0.0/src/main/native/ \
             -cp "${BISULFITEHIC_DIR}/target/bisulfitehic-0.38-jar-with-dependencies.jar:${BISULFITEHIC_DIR}/lib/jbwa.jar" \
             main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem \
             ${REFERENCE_GENOME} \
             ${OUTPUT_BAM} \
             ${FASTQ_R1} \
             ${FASTQ_R2} \
             -t 1 \
             -rgId ${SRR} \
             -rgSm sample_id \
             -outputMateDiffChr \
             -buffer 1000 \
             -nonDirectional \
             -pbat \
             -enzymeList ${ENZYME_LIST}

        # Check for errors in the Java execution
        if [ $? -ne 0 ]; then
            echo "Error: Bisulfite Hi-C mapping failed for ${SRR}."
            echo "${SRR}" >> $FAILED_LIST
            return 1
        fi
    else
        echo "Output BAM for ${SRR} already exists. Skipping mapping."
    fi

}

export -f process_srr

# Use GNU Parallel to download, convert, verify, and map SRA files in parallel
cat SRR_Acc_List.txt | parallel -j 8 process_srr

echo "All downloads, conversions, verifications, and mappings completed."


