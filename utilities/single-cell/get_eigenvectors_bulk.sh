#!/bin/bash

# Run the Python script and source the output to import the variables
eval $(python3 config_and_print.py)

# Ensure that the chromosomes variable is split into an array
IFS=',' read -r -a chromosomes <<< "${chromosomes}"

# Download the GM12878 hic file if it doesn't already exist
hic_GM12878_file="../../bin/softwarefiles/GSE63525_GM12878_30.hic"
if [ ! -f "$hic_GM12878_file" ]; then
    wget $hic_GM12878_url -O "$hic_GM12878_file"
fi

# Download the IMR90 hic file if it doesn't already exist
hic_IMR90_file="../../bin/softwarefiles/GSE63525_IMR90_30.hic"
if [ ! -f "$hic_IMR90_file" ]; then
    wget $hic_IMR90_url -O "$hic_IMR90_file"
fi

path_to_jar="../../bin/softwarefiles/juicer_tools_1.22.01.jar"

# Check and run the Java commands only if the files do not already exist
for res_label in ${resolutions[@]}; do
    resolution=$(echo $res_label | cut -d':' -f1)  # Extract just the numeric part of the resolution
    for chromosome in "${chromosomes[@]}"; do
        for data in ${data_types[@]}; do
            for prefix in GM12878 IMR90; do
                if [ "$prefix" == "GM12878" ]; then
                    path_to_hic="$hic_GM12878_file"
                else
                    path_to_hic="$hic_IMR90_file"
                fi

                eigenvector_directory="${output_directory}/eigenvector/"
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
done



