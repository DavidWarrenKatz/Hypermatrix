#!/bin/bash
# File: single_cell_pipeline.sh

##############################################################################
# Description: This script converts BAM files into Hi-C matrices and
# methylation matrices organized in folders for each chromosome.
# This pipeline also computes compartments and clustering and outputs
# images for comparison.
#
# In this phase, we are testing up to the make_hic_matrices.py script.
##############################################################################

set -euo pipefail

# Enable debugging if needed
# set -x

# Functions
check_path() {
    # Function to check if a file path is resolved correctly
    local script_name="$1"
    local path="$2"

    echo "[DEBUG]: Checking if the path for $script_name is valid..."
    if [ ! -f "$path" ]; then
        echo "[ERROR]: Failed to retrieve the path to $script_name (Path: $path)"
        exit 1
    else
        echo "[LOG]: Successfully resolved path for $script_name"
        echo " PATH: $path"
    fi
}

run_script() {
    # Function to execute a Python script and check its exit status
    local script_path="$1"
    local script_desc="$2"

    echo "[LOG]: Running $script_desc (Path: $script_path)"
    echo "[DEBUG]: Current HYPERMATRIX_CONFIG is $HYPERMATRIX_CONFIG"
    python "$script_path"

    local exit_code=$?
    echo "[DEBUG]: Exit code for $script_desc is $exit_code"

    if [ $exit_code -eq 0 ]; then
        echo "[LOG]: $script_desc completed successfully"
    else
        echo "[ERROR]: $script_desc failed to execute with exit code $exit_code"
        exit $exit_code
    fi
}

# Resolve HYPERMATRIX_DIR using the script's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYPERMATRIX_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"
echo "[DEBUG]: HYPERMATRIX_DIR resolved to $HYPERMATRIX_DIR"

# Set the configuration file path and export HYPERMATRIX_CONFIG
export HYPERMATRIX_CONFIG="${HYPERMATRIX_CONFIG:-$HYPERMATRIX_DIR/config.json}"
echo "[DEBUG]: HYPERMATRIX_CONFIG is set to $HYPERMATRIX_CONFIG"

# Construct paths to the scripts
EXPORT_CONFIG_PATH="$HYPERMATRIX_DIR/export_config.py"
MAKE_METHY_MATRICES_PATH="$HYPERMATRIX_DIR/utilities/single-cell/standard_pipeline/make_methy_matrices.py"
MAKE_HIC_MATRICES_PATH="$HYPERMATRIX_DIR/utilities/single-cell/standard_pipeline/make_hic_matrices.py"

# Check paths using the check_path function
check_path "export_config.py" "$EXPORT_CONFIG_PATH"
check_path "make_methy_matrices.py" "$MAKE_METHY_MATRICES_PATH"
check_path "make_hic_matrices.py" "$MAKE_HIC_MATRICES_PATH"

# Run the scripts and check exit status
run_script "$EXPORT_CONFIG_PATH" "Export configuration"
run_script "$MAKE_METHY_MATRICES_PATH" "Make methylation matrices"
run_script "$MAKE_HIC_MATRICES_PATH" "Make Hi-C matrices"

# Commented out other script executions for this phase
# echo "[LOG]: Skipping other scripts for this phase."

# Exit successfully
echo "[LOG]: single_cell_pipeline.sh completed successfully."
exit 0
