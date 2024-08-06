#!/bin/bash
# Description: The following script converts bam files into hic matrices organized in folders for each chromosome

# Run the Python script and source the output to import the variables
# each individual file still needs to import these variables as well
eval "$(python3 config_and_print.py)"

# Execute the filtering script
chmod +x filter_bam.sh
./filter_bam.sh

chmod +x filter_bw.sh
./filter_bw.sh

chmod +x update_filted_list.sh
./update_filted_list.sh

# Execute the make pairwise contact script
chmod +x make_pairwise_contact_format_from_bam.sh
./make_pairwise_contact_format_from_bam.sh

# Execute the filtering script
chmod +x preprocess_for_hic_cluster.sh
./preprocess_for_hic_cluster.sh

# Execute the make juicer format script
chmod +x make_juicer_short_format.sh
./make_juicer_short_format.sh

if [ "$cluster_compartments" = "True" ]; then
chmod +x compartment_calling.sh
./compartment_calling.sh
fi


