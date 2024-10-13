#!/bin/bash
# file: single-cell/preprocess/nomehic_preprocess.sh
# To-do 
# add shebang line 
# import parameters from config.json instead of config.py



# Get the directory where this script (single_cell_pipeline.sh) is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Import the parameters from config.py (relative to the script's directory)
eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

# Execute the filtering script (relative to the script's directory)
chmod +x $SCRIPT_DIR/filter_bam.sh
$SCRIPT_DIR/filter_bam.sh

# Generate the interval bed file
python $SCRIPT_DIR/generate_interval_bed.py

# Execute the filter bigwig files script
chmod +x $SCRIPT_DIR/filter_bw.sh
$SCRIPT_DIR/filter_bw.sh

# Execute the script that filters the filtered hic to contian high quality methy
#[TO DO: need to update this, right now hard-coded]
chmod +x $SCRIPT_DIR/update_filted_list.sh
$SCRIPT_DIR/update_filted_list.sh

# Execute the make pairwise contact script
chmod +x $SCRIPT_DIR/make_pairwise_contact_format_from_bam.sh
$SCRIPT_DIR/make_pairwise_contact_format_from_bam.sh

# Execute the preprocessing for hicluster script
chmod +x $SCRIPT_DIR/preprocess_for_hic_cluster.sh
$SCRIPT_DIR/preprocess_for_hic_cluster.sh

# Execute the make juicer format script
# This script forms the raw directory and the
#imputed directories
chmod +x $SCRIPT_DIR/make_juicer_short_format.sh
$SCRIPT_DIR/make_juicer_short_format.sh
