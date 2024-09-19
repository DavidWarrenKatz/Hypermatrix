#!/bin/bash

# Check if the 'bisulfitehic' environment exists
if conda env list | grep -q 'bisulfitehic'; then
    echo "Conda environment 'bisulfitehic' already exists. Skipping creation and setup."
else
    # Clone the repository from Bitbucket
    git clone https://bitbucket.org/dnaase/bisulfitehic.git

    cd bisulfitehic

    # Set the JAVA_HOME environment variable
    export JAVA_HOME=$(/usr/libexec/java_home)

    # Source the install.sh script
    source ./install.sh

    # Create a new conda environment named 'bisulfitehic' with Python 3.8
    conda create -n bisulfitehic python=3.8 -y

    # Activate the conda environment
    source activate bisulfitehic

    conda install -c bioconda bowtie2=2.4.2
    conda install -c bioconda bismark
    conda install -c bioconda bwa
    conda install -c bioconda samtools
    conda install -c bioconda picard
    pip install numpy
    pip install pysam

    echo "Setup complete. The 'bisulfitehic' environment is activated."
fi

# Directory where the mm9 genome should be downloaded
GENOME_DIR="../../projects/methyHic"
GENOME_FILE="$GENOME_DIR/mm9.fa"
FAI_FILE="$GENOME_DIR/mm9.fa.fai"
DICT_FILE="$GENOME_DIR/mm9.dict"

# Check if the mm9 genome has already been downloaded
if [ -f "$GENOME_FILE" ]; then
    echo "mm9 genome already exists. Skipping download."
else
    echo "Downloading mm9 genome..."
    wget -P "$GENOME_DIR" http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
    gunzip "$GENOME_DIR/mm9.fa.gz"
    echo "Download and extraction of mm9 genome complete."
fi

# Generate the .fai file using samtools if it doesn't exist
if [ ! -f "$FAI_FILE" ]; then
    echo "Generating .fai index file for mm9..."
    samtools faidx "$GENOME_FILE"
    echo "mm9.fa.fai generated."
else
    echo "mm9.fa.fai already exists. Skipping."
fi

# Generate the .dict file using Picard if it doesn't exist
if [ ! -f "$DICT_FILE" ]; then
    echo "Generating .dict file for mm9..."
    picard CreateSequenceDictionary REFERENCE="$GENOME_FILE" OUTPUT="$DICT_FILE"
    echo "mm9.dict generated."
else
    echo "mm9.dict already exists. Skipping."
fi


<<comment
# Check if the Bisulfite_Genome directory exists
if [ -d "$GENOME_DIR/Bisulfite_Genome" ]; then
    echo "Bisulfite Genome already prepared. Skipping bismark genome preparation."
else
    echo "Running bismark genome preparation..."
    bismark_genome_preparation "$GENOME_DIR"
    echo "Bismark genome preparation complete."
fi

# Index the bisulfite converted genomes using bwa
CT_CONVERSION_FA="$GENOME_DIR/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
GA_CONVERSION_FA="$GENOME_DIR/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"

if [ -f "$CT_CONVERSION_FA.bwt" ] && [ -f "$GA_CONVERSION_FA.bwt" ]; then
    echo "BWA index for bisulfite converted genomes already exists. Skipping indexing."
else
    echo "Indexing bisulfite converted genomes using bwa..."
    bwa index "$CT_CONVERSION_FA"
    bwa index "$GA_CONVERSION_FA"
    echo "BWA indexing complete."
fi

# Confirm expected files are present
echo "Checking final file structure in the reference genome folder:"
if [ -f "$GENOME_FILE" ]; then
    echo "mm9.fa exists."
else
    echo "mm9.fa is missing!"
fi

if [ -f "$GENOME_FILE.fai" ]; then
    echo "mm9.fa.fai exists."
else
    echo "mm9.fa.fai is missing!"
fi

if [ -f "$GENOME_DIR/mm9.dict" ]; then
    echo "mm9.dict exists."
else
    echo "mm9.dict is missing!"
fi

if [ -d "$GENOME_DIR/Bisulfite_Genome/CT_conversion" ] && [ -f "$CT_CONVERSION_FA" ]; then
    echo "CT_conversion genome and BWA index files exist."
else
    echo "CT_conversion genome or BWA index files are missing!"
fi

if [ -d "$GENOME_DIR/Bisulfite_Genome/GA_conversion" ] && [ -f "$GA_CONVERSION_FA" ]; then
    echo "GA_conversion genome and BWA index files exist."
else
    echo "GA_conversion genome or BWA index files are missing!"
fi

comment
