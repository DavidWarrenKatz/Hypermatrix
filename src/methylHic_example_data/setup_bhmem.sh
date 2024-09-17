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

    # Create a new conda environment named 'bisulfitehic' with Python 3.6
    conda create -n bisulfitehic python=3.6 -y

    # Activate the conda environment
    source activate bisulfitehic

    # Install Bismark
    conda install -c bioconda bismark

    # Install the required Python packages
    pip install numpy
    pip install pysam

    echo "Setup complete. The 'bisulfitehic' environment is activated."
fi

# Directory where the mm9 genome should be downloaded
GENOME_DIR="../../projects/methyHic"
GENOME_FILE="$GENOME_DIR/mm9.fa"

# Check if the mm9 genome has already been downloaded
if [ -f "$GENOME_FILE" ]; then
    echo "mm9 genome already exists. Skipping download."
else
    echo "Downloading mm9 genome..."
    wget -P "$GENOME_DIR" http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
    gunzip "$GENOME_DIR/mm9.fa.gz"
    echo "Download and extraction of mm9 genome complete."
fi

