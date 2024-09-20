#!/bin/bash

# Load necessary modules
module load sra-toolkit/2.10.9
module load FastQC/0.11.9

# Create necessary directories
mkdir -p ../../projects/methyHic
mkdir -p ../../projects/methyHic/fastqc_reports
mkdir -p ../../projects/methyHic/output_bam  # Directory for output BAM files

# Files to store success and failure lists
SUCCESS_LIST="../../projects/methyHic/successful_fastqs.txt"
FAILED_LIST="../../projects/methyHic/failed_fastqs.txt"

# Initialize or clear the success and failure lists
echo -n > $SUCCESS_LIST
echo -n > $FAILED_LIST

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

# Function to download, convert, verify each SRR, and run bisulfite Hi-C mapping
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

    # Check if FastQC report already exists
    FASTQC_REPORT="../../projects/methyHic/fastqc_reports/${SRR}_1_fastqc.html"
    if [ ! -f "${FASTQC_REPORT}" ]; then
        # Run FastQC to check the integrity of the FASTQ files
        echo "Running FastQC on ${SRR} FASTQ files..."
        fastqc ../../projects/methyHic/${SRR}_*.fastq.gz -o ../../projects/methyHic/fastqc_reports/
        if [ $? -ne 0 ]; then
            echo "Error: FastQC failed for ${SRR}."
            echo "${SRR}" >> $FAILED_LIST
            return 1
        fi
    else
        echo "FastQC report for ${SRR} already exists. Skipping FastQC."
    fi

    # Check if output BAM file already exists
    OUTPUT_BAM="../../projects/methyHic/output_bam/${SRR}_output.bam"
    if [ ! -f "${OUTPUT_BAM}" ]; then
        # Run the bisulfite Hi-C mapping using Java
        echo "Running bisulfite Hi-C mapping for ${SRR}..."

        # Set variables for the paths to make the command more readable
        BISULFITEHIC_DIR="./bisulfitehic"
        REFERENCE_GENOME="../../projects/methyHic/mm9.fa"
        FASTQ_R1="../../projects/methyHic/${SRR}_1.fastq.gz"
        FASTQ_R2="../../projects/methyHic/${SRR}_2.fastq.gz"
        ENZYME_LIST="../../projects/methyHic/mm9_DpnII.bed"

        # Java command for bisulfite Hi-C mapping
        java -Xmx15G \
             -Djava.library.path=${BISULFITEHIC_DIR}/jbwa/src/main/native/ \
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

    # If everything was successful, add to success list
    echo "${SRR}" >> $SUCCESS_LIST

    echo "[$(date)] Finished processing ${SRR} with PID $$"
}

export -f process_srr

# Use GNU Parallel to download, convert, verify, and map SRA files in parallel
cat SRR_Acc_List.txt | parallel -j 8 process_srr

echo "All downloads, conversions, verifications, and mappings completed."

# Summary of results
echo "Summary of FASTQ processing:"
echo "Successfully processed FASTQs: $(wc -l < $SUCCESS_LIST)"
echo "Failed FASTQs: $(wc -l < $FAILED_LIST)"

