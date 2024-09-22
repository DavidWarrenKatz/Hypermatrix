#!/bin/bash
# File: test_single_cell_pipeline.sh

##############################################################################
# Description: This script tests the single_cell_pipeline.sh script.
# It sets up the environment, runs the pipeline, and checks outputs.
#
# In this phase, we are testing up to the make_hic_matrices.py script.
##############################################################################

set -euo pipefail  # Keep 'set -u' for strict variable checking

# Enable debugging if needed
# set -x

# Temporarily disable 'set -u' before initializing Conda
set +u
echo "[DEBUG]: Initializing Conda for the current shell..."
eval "$(conda shell.bash hook)"

# Activate the hypermatrix Conda environment
echo "[DEBUG]: Activating hypermatrix Conda environment..."
conda activate hypermatrix

# Re-enable 'set -u' after activation
set -u

# Set variables
TEST_DATA_DIR="$HOME/test_data"
TEST_OUTPUT_DIR="$TEST_DATA_DIR/output"
HYPERMATRIX_DIR="$HOME/Hypermatrix/hypermatrix"

echo "[DEBUG]: TEST_DATA_DIR is set to $TEST_DATA_DIR"
echo "[DEBUG]: TEST_OUTPUT_DIR is set to $TEST_OUTPUT_DIR"
echo "[DEBUG]: HYPERMATRIX_DIR is set to $HYPERMATRIX_DIR"

# Ensure output directory exists
mkdir -p "$TEST_OUTPUT_DIR"

# Set the configuration file path and export HYPERMATRIX_CONFIG
export HYPERMATRIX_CONFIG="$TEST_OUTPUT_DIR/config.json"
echo "[DEBUG]: HYPERMATRIX_CONFIG is set to $HYPERMATRIX_CONFIG"

# Copy or create a minimal config.json if necessary
if [ ! -f "$HYPERMATRIX_CONFIG" ]; then
    echo "[DEBUG]: config.json not found at $HYPERMATRIX_CONFIG. Copying default config."
    cp "$HYPERMATRIX_DIR/config.json" "$HYPERMATRIX_CONFIG"
else
    echo "[DEBUG]: config.json already exists at $HYPERMATRIX_CONFIG"
fi

# Update the config.json with test-specific paths
echo "[DEBUG]: Updating config.json with test-specific paths..."
python3 - << EOF
import json
config_path = "$HYPERMATRIX_CONFIG"
with open(config_path, 'r') as f:
    config = json.load(f)

config['bam_directory'] = "$TEST_DATA_DIR"
config['methy_directory'] = "$TEST_DATA_DIR"
config['output_directory'] = "$TEST_OUTPUT_DIR"
config['genome_fa'] = "$TEST_DATA_DIR/hg38.fa.gz"
config['chrom_file'] = "$TEST_DATA_DIR/hg38.chrom.sizes"
config['filtered_list'] = "$TEST_OUTPUT_DIR/filtered_bam_list.txt"
config['software_directory'] = "$HYPERMATRIX_DIR/genomes"
config['fragments_file'] = "$HYPERMATRIX_DIR/genomes/hg38_DpnII.txt"
config['resolutions'] = ["1000000:1Mb"]  # Adjusted resolution
config['chromosomes'] = [f"chr{i}" for i in range(1, 23)]  # Include 'chr' prefix
config['methy_data_file'] = "GSM1386019_oocyte_mc_CG_plus.bed.gz"
config['hic_cell_type1_url'] = ""  # Update with actual URL if needed
config['hic_cell_type2_url'] = ""  # Update with actual URL if needed

with open(config_path, 'w') as f:
    json.dump(config, f, indent=4)
EOF

echo "[DEBUG]: Updated config.json content:"
cat "$HYPERMATRIX_CONFIG"

# Export HYPERMATRIX_DIR for use in the pipeline script
export HYPERMATRIX_DIR

# Run the single_cell_pipeline.sh script
PIPELINE_SCRIPT="$HYPERMATRIX_DIR/utilities/single-cell/standard_pipeline/single_cell_pipeline.sh"
echo "[DEBUG]: Running single_cell_pipeline.sh script at $PIPELINE_SCRIPT"
"$PIPELINE_SCRIPT"

# For this phase, we are testing up to make_hic_matrices.py
# We will check if the Hi-C matrices have been generated correctly

echo "[DEBUG]: Checking generated Hi-C matrices..."

# Define the expected output directory for Hi-C matrices
HIC_MATRICES_DIR="$TEST_OUTPUT_DIR/hic_matrices"

# Check if the directory exists and contains expected files
if [ -d "$HIC_MATRICES_DIR" ] && [ "$(ls -A "$HIC_MATRICES_DIR")" ]; then
    echo "[PASS] Hi-C matrices generated successfully."
else
    echo "[FAIL] Hi-C matrices not generated."
    exit 1
fi

echo "[LOG]: Test of make_hic_matrices.py completed successfully."
