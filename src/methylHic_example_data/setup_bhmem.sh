#!/bin/bash
# file: setup_bhmem.sh



# Load conda functions into the shell
eval "$(conda shell.bash hook)"

# Function to check if a package is installed in the environment
function check_and_install_package() {
    # to-do add channel arg $2
    # alter to use bioconda channel
    package_name=$1
    if ! conda list | grep -q "^$package_name"; then
        echo "$package_name not found. Installing..."
        conda install -y -c bioconda -c conda-forge $package_name || { echo "Failed to install $package_name"; exit 1; }
    else
        echo "$package_name is already installed."
    fi
}

# Check if the 'bisulfitehic' environment exists
if conda env list | grep -q 'bisulfitehic'; then
    echo "Conda environment 'bisulfitehic' already exists."
    source $(conda info --base)/etc/profile.d/conda.sh
     
    
    # Check and install necessary packages if not already installed
    check_and_install_package openjdk
    check_and_install_package gcc
    check_and_install_package bismark
    check_and_install_package picard
    check_and_install_package bwa
    
    # Check and install Python packages via pip
    pip show numpy || pip install numpy
    pip show pysam || pip install pysam


    # ensure bitbucket exist if not download it 
    # Download the directory
    #wget https://bitbucket.org/dnaase/bisulfitehic/get/04680506dd40.zip || { echo "Failed to download zip file"; exit 1; }
    #zip_file="04680506dd40.zip"
    output_folder="bisulfitehic"


    # Unzip the file into the custom output folder
    # unzip "$zip_file" -d "$output_folder" || { echo "Failed to unzip file"; exit 1; }
    # rm "$zip_file" # Clean up zip file

 

else
    # Create the conda environment if it doesn't exist
    echo "Creating conda environment 'bisulfitehic'..."
    conda create -n bisulfitehic python=3.6 -y || { echo "Conda environment creation failed!"; exit 1; }
    
    # Activate the environment
    conda activate bisulfitehic
    
    # Install packages
    conda install -y -c conda-forge openjdk gcc
    conda install -y -c bioconda bismark picard bwa
    
    # Install Python packages via pip
    pip install numpy pysam
fi

# Additional setup commands if needed (e.g., set JAVA_HOME)
java_path=$(which java)
if [ -z "$java_path" ]; then
    echo "Java not found in the current Conda environment."
    exit 1
fi

java_home=$(dirname $(dirname "$java_path"))
if [ -d "$java_home/include" ]; then
    export JAVA_HOME="$java_home"
    echo "JAVA_HOME is set to: $JAVA_HOME"
else
    echo "The Java installation found is not a JDK. Please install a JDK."
    exit 1
fi

pwd
cd bisulfitehic/dnaase-bisulfitehic-04680506dd40 || { echo "Failed to change directory"; exit 1; }
pwd

echo "running install.sh"
#  Source the install.sh script to install bisulfitehic software
source ./install.sh || { echo "Failed to install bisulfitehic"; exit 1; }


# Genome directory
GENOME_DIR="${1:-"../../../../projects/methyHic"}"
mkdir -p "$GENOME_DIR"

echo "[LOG]: showing working directory $GENOME_DIR"

GENOME_FILE="$GENOME_DIR/mm9.fa"
FAI_FILE="$GENOME_FILE.fai"
DICT_FILE="$GENOME_DIR/mm9.dict"

# mm9 genome has already been downloaded?
if [ -f "$GENOME_FILE" ]; then
    echo "mm9 genome already exists. Skipping download."
else
    echo "Downloading mm9 genome..."
    wget -P "$GENOME_DIR" http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz || { echo "Failed to download mm9 genome"; exit 1; }
    gunzip "$GENOME_DIR/mm9.fa.gz" || { echo "Failed to extract mm9 genome"; exit 1; }
    echo "Download and extraction of mm9 genome complete."
fi

# Generate the .fai file if it doesn't exist
if [ ! -f "$FAI_FILE" ]; then
    echo "Generating .fai index file for mm9..."
    samtools faidx "$GENOME_FILE" || { echo "Failed to generate .fai file"; exit 1; }
else
    echo ".fai file already exists."
fi

# Generate the .dict file if it doesn't exist
if [ ! -f "$DICT_FILE" ]; then
    echo "Generating .dict file for mm9..."
    picard CreateSequenceDictionary REFERENCE="$GENOME_FILE" OUTPUT="$DICT_FILE" || { echo "Failed to generate .dict file"; exit 1; }
else
    echo ".dict file already exists."
fi

# Prepare Bisulfite Genome
if [ -d "$GENOME_DIR/Bisulfite_Genome" ]; then
    echo "Bisulfite Genome already prepared. Skipping bismark genome preparation."
else
    echo "Running bismark genome preparation..."
    bismark_genome_preparation --parallel 64 "$GENOME_DIR" || { echo "Failed to prepare bisulfite genome"; exit 1; }
fi

# Index bisulfite genomes using bwa
CT_CONVERSION_FA="$GENOME_DIR/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
GA_CONVERSION_FA="$GENOME_DIR/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"

if [ -f "$CT_CONVERSION_FA.bwt" ] && [ -f "$GA_CONVERSION_FA.bwt" ]; then
    echo "BWA index for bisulfite genomes already exists."
else
    echo "Indexing bisulfite genomes with bwa..."
    bwa index "$CT_CONVERSION_FA" || { echo "Failed to index CT genome"; exit 1; }
    bwa index "$GA_CONVERSION_FA" || { echo "Failed to index GA genome"; exit 1; }
    echo "BWA indexing complete."
fi

# Confirm expected files are present
echo "Final file structure in the genome folder:"
[ -f "$GENOME_FILE" ] && echo "mm9.fa exists." || echo "mm9.fa is missing!"
[ -f "$FAI_FILE" ] && echo ".fai file exists." || echo ".fai file is missing!"
[ -f "$DICT_FILE" ] && echo ".dict file exists." || echo ".dict file is missing!"
[ -d "$GENOME_DIR/Bisulfite_Genome/CT_conversion" ] && [ -f "$CT_CONVERSION_FA" ] && echo "CT genome files exist." || echo "CT genome files are missing!"
[ -d "$GENOME_DIR/Bisulfite_Genome/GA_conversion" ] && [ -f "$GA_CONVERSION_FA" ] && echo "GA genome files exist." || echo "GA genome files are missing!"

echo "Setup and indexing complete."
