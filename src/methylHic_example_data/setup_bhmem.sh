#!/bin/bash

# Check if the 'bisulfitehic' environment exists
if conda env list | grep -q 'bisulfitehic'; then
    echo "Conda environment 'bisulfitehic' already exists. Skipping creation and setup."
    source activate bisulfitehic
else
    #Create the conda bisulfitehic environment according to yml 
    conda env create -f bisulfitehic_environment.yml

    # Activate the conda environment
    source activate bisulfitehic

    # Download the directory
    wget https://bitbucket.org/dnaase/bisulfitehic/get/04680506dd40.zip
    zip_file="04680506dd40.zip"
    output_folder=bisulfitehic  

    # Create the output folder if it doesn't exist
    mkdir -p "$output_folder"

    # Unzip the file into the custom output folder
    unzip "$zip_file" -d "$output_folder"

    # Dynamically set the JAVA_HOME environment variable
    # Find the Java executable within the Conda environment
    java_path=$(which java)

    # Check if Java is installed
    if [ -z "$java_path" ]; then
        echo "Java not found in the current Conda environment."
        exit 1
    fi

    # Get the directory of the Java installation
    java_home=$(dirname $(dirname "$java_path"))

    # Set JAVA_HOME and export it
    export JAVA_HOME="$java_home"

    # Verify JAVA_HOME is set correctly
    if [ -z "$JAVA_HOME" ]; then
        echo "Failed to set JAVA_HOME."
    else
        echo "JAVA_HOME is set to: $JAVA_HOME"
    fi

    cd bisulfitehic/dnaase-bisulfitehic-04680506dd40

    # Source the install.sh script to install bisulfitehic software into enviroment
    source ./install.sh

    # install the remianing packages
    conda install -c bioconda bowtie2=2.4.2
    conda install -c bioconda bismark
    conda install -c bioconda bwa
    conda install -c bioconda samtools
    conda install -c bioconda picard
    pip install numpy
    pip install pysam

    echo "Setup complete. The 'bisulfitehic' environment is activated."

    cd ..
    cd ..
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


