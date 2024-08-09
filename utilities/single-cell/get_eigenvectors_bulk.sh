#!/bin/bash

# Run the Python script and source the output to import the variables
eval $(python3 config_and_print.py)

# Download the GM12878 hic file if it doesn't already exist
hic_GM12878_file="../../bin/softwarefiles/GSE63525_GM12878.hic"
if [ ! -f "$dark_regions_file" ]; then
    wget $hic_GM12878_url -O "$hic_GM12878_file"
fi

# Download the IMR90 hic file if it doesn't already exist
hic_IMR90_file="../../bin/softwarefiles/GSE63525_IMR90.hic"
if [ ! -f "$dark_regions_file" ]; then
    wget $hic_IMR90_url -O "$hic_IMR90_file"
fi

path_to_jar="../projects/softwarefiles/juicer_tools_1.22.01.jar"

# Check and run the Java commands only if the files do not already exist
for resolution in ${resolutions[@]}; do
    for chromosome in ${chromosomes[@]}; do
        for data in ${data_types[@]}; do
            for prefix in ['GM12878', 'IMR90']
                path_to_new_hic="../../bin/softwarefiles/GSE63525_${prefix}IMR90.hic"
                eigenvector_directory="${output_directory}eigenvector/"
                mkdir -p $eigenvector_directory
                path_to_new_eigenvector="${eigenvector_directory}res${resolution}_ch${chromosome}_${data}_${prefix}_KR_eigenvector.txt"

                if [[ ! -f "$path_to_new_eigenvector" ]]; then
                    eigenvector_command="java -jar ${path_to_jar} eigenvector -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_eigenvector}"
                    echo "Running: $eigenvector_command"
                    $eigenvector_command
                else
                    echo "Eigenvector file already exists: $path_to_new_eigenvector"
                fi
        done
    done
done	
