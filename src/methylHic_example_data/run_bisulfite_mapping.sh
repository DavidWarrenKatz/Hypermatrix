#!/bin/bash

# Create necessary directories
mkdir -p ../../projects/methyHic/output_bam  # Directory for output BAM files

# Function to run bisulfite Hi-C mapping, process BAM with samtools, Bis-SNP, and Juicer
run_mapping() {
    SRR=$1
    echo "[$(date)] Running bisulfite Hi-C mapping for ${SRR}"

    OUTPUT_BAM="../../projects/methyHic/output_bam/${SRR}_output.bam"

    # Check if output BAM file already exists
    if [ ! -f "${OUTPUT_BAM}" ]; then
        # Set variables for the paths to make the command more readable
        BISULFITEHIC_DIR="./bisulfitehic/dnaase-bisulfitehic-04680506dd40"
        REFERENCE_GENOME="../../projects/methyHic/mm9.fa"
        FASTQ_R1="../../projects/methyHic/${SRR}_1.fastq.gz"
        FASTQ_R2="../../projects/methyHic/${SRR}_2.fastq.gz"
        ENZYME_LIST="../../projects/methyHic/mm9_DpnII.bed"

        # Run the bisulfite Hi-C mapping using Java
        java -Xmx15G \
             -Djava.library.path=${BISULFITEHIC_DIR}/jbwa/jbwa-1.0.0/src/main/native \
             -cp "${BISULFITEHIC_DIR}/target/bisulfitehic-0.38-jar-with-dependencies.jar:${BISULFITEHIC_DIR}/jbwa/jbwa-1.0.0/jbwa.jar" \
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
            return 1
        fi

        # Samtools processing steps after bisulfite Hi-C mapping
        echo "[$(date)] Running samtools processing for ${SRR}..."

        # Define input and output BAM files
        SORTED_BAM="../../projects/methyHic/output_bam/${SRR}_sorted.bam"
        CALMD_BAM="../../projects/methyHic/output_bam/${SRR}_calmd.bam"

        # Samtools commands for sorting, fixing, marking duplicates, and recalibrating base qualities
        samtools sort --threads 5 -T ${SRR}.name -n ${OUTPUT_BAM} | \
        samtools fixmate -m --threads 5 - - | \
        samtools sort --threads 5 -T ${SRR}.cor - | \
        samtools markdup -T ${SRR}.mdups --threads 5 - - | \
        samtools calmd --threads 5 -b - ${REFERENCE_GENOME} 2>/dev/null > ${CALMD_BAM}

        # Index the final BAM file
        samtools index -@ 10 ${CALMD_BAM}

        echo "Samtools processing for ${SRR} completed and output BAM is indexed."

        # Run Bis-SNP after samtools processing
        echo "[$(date)] Running Bis-SNP for ${SRR}..."

        # Define paths to Bis-SNP, reference genome, and dbSNP
        BISSNP_DIR="/your/path/to/Bis-SNP"
        REFERENCE_GENOME="/your/path/to/reference_genome/human_g1k_v37.fa"
        DBSNP_VCF="/your/path/to/dbsnp/dbsnp_138.b37.vcf"
        BISSNP_JAR="${BISSNP_DIR}/BisSNP-0.90.jar"

        # Run Bis-SNP using Perl script
        perl ${BISSNP_DIR}/bissnp_easy_usage.pl \
             --use_bad_mates \
             --mmq 30 \
             --nt 30 \
             --mem 5 \
             ${BISSNP_JAR} \
             ${CALMD_BAM} \
             ${REFERENCE_GENOME} \
             ${DBSNP_VCF}

        echo "Bis-SNP processing for ${SRR} completed."

        # Run sam2juicer and Juicer steps
        echo "[$(date)] Running sam2juicer and Juicer processing for ${SRR}..."

        # Define paths to sam2juicer and Juicer tools
        SAM2JUICER_PY="/your/path/to/bisulfitehic/src/python/sam2juicer.py"
        JUICER_RESTRICTION_SITES="/your/path/to/juicer/restriction_sites/hg19_DpnII.txt"
        JUICER_JAR="/your/path/to/juicer/scripts/common/juicer_tools.jar"

        # Run sam2juicer
        python ${SAM2JUICER_PY} -s ${CALMD_BAM} -f ${JUICER_RESTRICTION_SITES} | gzip -fc > ${SRR}.hic.txt.gz

        # Run Juicer tools to generate the .hic file
        java -Xmx20G -jar ${JUICER_JAR} pre -d -q 30 -f ${JUICER_RESTRICTION_SITES} ${SRR}.hic.txt.gz ${SRR}.mapQ30.hic hg19

        echo "sam2juicer and Juicer processing for ${SRR} completed."
    else
        echo "Output BAM for ${SRR} already exists. Skipping mapping, samtools, Bis-SNP, and Juicer processing."
    fi
}

export -f run_mapping

# Use GNU Parallel to run bisulfite Hi-C mapping, samtools processing, Bis-SNP, and Juicer for each SRR in parallel
cat SRR_Acc_List.txt | parallel -j 8 run_mapping

echo "All bisulfite Hi-C mappings, samtools processing, Bis-SNP, and Juicer analysis completed."

