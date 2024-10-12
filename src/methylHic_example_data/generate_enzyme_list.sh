#!/bin/bash
# file: generate_enzyme_list.sh

# Enable debugging for the script
# set -x  # Turn on tracing to see each command executed

# Define the directory to install HiC-Pro and the output file for the enzyme list
HICPRO_DIR="../../softwarefiles/HiC-Pro"  # Path for HiC-Pro installation
UTIL_DIR="$HICPRO_DIR/bin/utils"
MM9_GENOME="../../../../projects/methyHic/mm9.fa"
OUTPUT_ENZYME_LIST="../../../../projects/methyHic/mm9_DpnII.bed"

# Debug: Trace the directories and paths
echo "HICPRO_DIR is set to: $HICPRO_DIR"
echo "UTIL_DIR is set to: $UTIL_DIR"
echo "MM9_GENOME is set to: $MM9_GENOME"
echo "OUTPUT_ENZYME_LIST is set to: $OUTPUT_ENZYME_LIST"

# Check if HiC-Pro is installed
if [ ! -d "$HICPRO_DIR" ]; then
    echo "HiC-Pro not found. Installing HiC-Pro..."

    # Clone the HiC-Pro repository
    git clone https://github.com/nservant/HiC-Pro.git $HICPRO_DIR
    
    # Debug: Check if the directory was created successfully
    if [ -d "$HICPRO_DIR" ]; then
        echo "HiC-Pro directory created successfully."
    else
        echo "Error: HiC-Pro directory was not created. Check permissions and paths."
        exit 1
    fi

    # Navigate to the HiC-Pro directory and install dependencies
    cd $HICPRO_DIR || { echo "Failed to change directory to $HICPRO_DIR"; exit 1; }
    make configure || { echo "Failed to configure HiC-Pro"; exit 1; }
    make install || { echo "Failed to install HiC-Pro"; exit 1; }
    
    echo "HiC-Pro installed successfully."
else
    echo "HiC-Pro already installed at $HICPRO_DIR"
fi

# Debug: Confirm the HiC-Pro directory and utility paths
if [ ! -d "$UTIL_DIR" ]; then
    echo "Error: Utilities directory $UTIL_DIR does not exist!"
    exit 1
else
    echo "Utilities directory exists: $UTIL_DIR"
fi

# Generate the DpnII enzyme list for mm9
if [ ! -f "$OUTPUT_ENZYME_LIST" ]; then
    echo "Generating DpnII enzyme list for mm9..."
    
    # Run the digest_genome.py script to generate the enzyme list
    python $UTIL_DIR/digest_genome.py -r dpnii -o $OUTPUT_ENZYME_LIST $MM9_GENOME || { echo "Failed to generate enzyme list"; exit 1; }
    
    echo "Enzyme list generated: $OUTPUT_ENZYME_LIST"
else
    echo "Enzyme list already exists: $OUTPUT_ENZYME_LIST"
fi

# Final confirmation
if [ -f "$OUTPUT_ENZYME_LIST" ]; then
    echo "Enzyme list generated successfully and located at: $OUTPUT_ENZYME_LIST"
else
    echo "Error: Enzyme list not found at: $OUTPUT_ENZYME_LIST"
    exit 1
fi

echo "All done!"

# Disable tracing after execution
set +x
