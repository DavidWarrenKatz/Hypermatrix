#!/bin/bash

# Define the directory to install HiC-Pro and the output file for the enzyme list
HICPRO_DIR="../../softwarefiles/HiC-Pro"  # Correct the path here
UTIL_DIR="$HICPRO_DIR/bin/utils"
MM9_GENOME="../../projects/methyHic/mm9.fa"
OUTPUT_ENZYME_LIST="../../projects/methyHic/mm9_DpnII.bed"

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
    
    # Run the digest_genome.py script to generate the enzyme list
    python $UTIL_DIR/digest_genome.py -r dpnii -o $OUTPUT_ENZYME_LIST $MM9_GENOME
    
    echo "Enzyme list generated: $OUTPUT_ENZYME_LIST"
else
    echo "Enzyme list already exists: $OUTPUT_ENZYME_LIST"
fi

echo "All done!"

