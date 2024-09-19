#!/bin/bash

# Define the directory to install HiC-Pro and the output file for the enzyme list
HICPRO_DIR="$../softwarefiles/HiC-Pro"
# Check if the directory exists, if not, create it
if [ ! -d "$HICPRO_DIR" ]; then
  mkdir -p "$HICPRO_DIR"
  echo "Directory $HICPRO_DIR created."
else
  echo "Directory $HICPRO_DIR already exists."
fi

MM9_GENOME="../../projects/methyHic/ m9.fa"
OUTPUT_ENZYME_LIST="mm9_DpnII.bed"

# Check if HiC-Pro is installed
if [ ! -d "$HICPRO_DIR" ]; then
    echo "HiC-Pro not found. Installing HiC-Pro..."

    # Clone the HiC-Pro repository
    git clone https://github.com/nservant/HiC-Pro.git $HICPRO_DIR
    
    # Navigate to the HiC-Pro directory and install dependencies
    cd $HICPRO_DIR
    make configure
    make install
    
    echo "HiC-Pro installed successfully."
else
    echo "HiC-Pro already installed."
fi

# Generate the DpnII enzyme list for mm9
if [ ! -f "$OUTPUT_ENZYME_LIST" ]; then
    echo "Generating DpnII enzyme list for mm9..."
    
    # Navigate to the HiC-Pro bin directory to run digest_genome.py
    cd $HICPRO_DIR/scripts
    
    # Run the digest_genome.py script to generate the enzyme list
    python digest_genome.py -r DpnII -o $OUTPUT_ENZYME_LIST -g $MM9_GENOME
    
    echo "Enzyme list generated: $OUTPUT_ENZYME_LIST"
else
    echo "Enzyme list already exists: $OUTPUT_ENZYME_LIST"
fi

echo "All done!"

