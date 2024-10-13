#!/bin/bash
# file: run_bisulfite_mapping.sh


# To-do
# see line 51 BISULFITEHIC_DIR does not exist. clone bit bucket

eval "$(conda shell.bash hook)"
conda activate bisulfitehic

# Define the path for genome download and extraction
GENOME_DIR="../src"
GENOME_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
GENOME_FA="${GENOME_DIR}/human_g1k_v37.fa"

# Download and unzip the genome if not already present
if [ ! -f "${GENOME_FA}" ]; then
    echo "[$(date)] Downloading the genome..."
    wget -P ${GENOME_DIR} ${GENOME_URL}
    gunzip ${GENOME_DIR}/human_g1k_v37.fasta.gz
    mv ${GENOME_DIR}/human_g1k_v37.fasta ${GENOME_FA}
    echo "[$(date)] Genome download and extraction completed."
else
    echo "[$(date)] Genome already exists. Skipping download."
fi

# Download dbSNP VCF file and unzip if necessary
DBSNP_DIR="../dbsnp"
mkdir -p ${DBSNP_DIR}
DBSNP_URL="https://data.broadinstitute.org/snowman/hg19/variant_calling/vqsr_resources/Exome/v2/dbsnp_138.b37.vcf.gz"
DBSNP_VCF="${DBSNP_DIR}/dbsnp_138.b37.vcf"

if [ ! -f "${DBSNP_VCF}" ]; then
    echo "[$(date)] Downloading dbSNP file..."
    wget -c ${DBSNP_URL} -P ${DBSNP_DIR}
    gunzip ${DBSNP_DIR}/dbsnp_138.b37.vcf.gz
    echo "[$(date)] dbSNP file download and extraction completed."
else
    echo "[$(date)] dbSNP file already exists. Skipping download."
fi

mkdir -p ../../../../projects/methyHic/output_bam

# Function to run bisulfite Hi-C mapping, process BAM with samtools, Bis-SNP, and Juicer
run_mapping() {
    SRR=$1

    OUTPUT_BAM="../../../../projects/methyHic/output_bam/${SRR}_output.bam"
    
    if [ ! -f "${OUTPUT_BAM}" ]; then
        BISULFITEHIC_DIR="./bisulfitehic/dnaase-bisulfitehic-04680506dd40"
        REFERENCE_GENOME="${GENOME_FA}"
        FASTQ_R1="../../../../projects/methyHic/${SRR}_1.fastq.gz"
        FASTQ_R2="../../../../projects/methyHic/${SRR}_2.fastq.gz"
        ENZYME_LIST="../../../../projects/methyHic/mm9_DpnII.bed"

        
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

        if [ $? -ne 0 ]; then
            echo "Error: Bisulfite Hi-C mapping failed for ${SRR}."
            return 1
        fi

        SORTED_BAM="../../../../projects/methyHic/output_bam/${SRR}_sorted.bam"
        CALMD_BAM="../../../../projects/methyHic/output_bam/${SRR}_calmd.bam"

        # Samtools commands for sorting, fixing, marking duplicates, and recalibrating base qualities
        samtools sort --threads 64 -T ${SRR}.name -n ${OUTPUT_BAM} | \
        samtools fixmate -m --threads 64 - - | \
        samtools sort --threads 64 -T ${SRR}.cor - | \
        samtools markdup -T ${SRR}.mdups --threads 5 - - | \
        samtools calmd --threads 64 -b - ${REFERENCE_GENOME} 2>/dev/null > ${CALMD_BAM}

        # Index the final BAM file
        samtools index -@ 64 ${CALMD_BAM}

        # Run Bis-SNP after samtools processing
        BISSNP_DIR="."
        BISSNP_JAR="${BISSNP_DIR}/BisSNP-0.90.jar"

        perl ${BISSNP_DIR}/bissnp_easy_usage.pl \
             --use_bad_mates \
             --mmq 30 \
             --nt 30 \
             --mem 5 \
             ${BISSNP_JAR} \
             ${CALMD_BAM} \
             ${REFERENCE_GENOME} \
             ${DBSNP_VCF}

        # Run sam2juicer and Juicer steps
        SAM2JUICER_PY="../sam2juicer.py"
        JUICER_RESTRICTION_SITES="../hg19_DpnII.txt"
        JUICER_JAR="../juicer_tools_1.22.01.jar"

        python ${SAM2JUICER_PY} -s ${CALMD_BAM} -f ${JUICER_RESTRICTION_SITES} | gzip -fc > ${SRR}.hic.txt.gz

        java -Xmx20G -jar ${JUICER_JAR} pre -d -q 30 -f ${JUICER_RESTRICTION_SITES} ${SRR}.hic.txt.gz ${SRR}.mapQ30.hic hg19
    else
        echo "Output BAM for ${SRR} already exists. Skipping mapping, samtools, Bis-SNP, and Juicer processing."
    fi
}

export -f run_mapping

# Use GNU Parallel to run bisulfite Hi-C mapping, samtools processing, Bis-SNP, and Juicer for each SRR in parallel
cat "./SRR_Acc_List.txt" | parallel -j 64 run_mapping

echo "All bisulfite Hi-C mappings, samtools processing, Bis-SNP, and Juicer analysis completed."
