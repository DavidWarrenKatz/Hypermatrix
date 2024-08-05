#!/bin/bash
# Description: The following script pulls bulk hi-c matrices for further processing downstream

# Variables
# Run the Python script and source the output to import the variables
eval $(python3 config_and_print.py)

# Create output directory if it doesn't exist
mkdir -p $output_directory

./filter_bam.sh


