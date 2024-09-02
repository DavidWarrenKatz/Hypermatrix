#!/bin/bash

# Run the Python script and source the output to import the variables
eval "$(python3 config_and_print.py)"

# Ensure that the chromosomes variable is treated as an array
#chromosomes=(${chromosomes[@]})

# Ensure that the chromosomes variable is split into an array
IFS=',' read -r -a chromosomes <<< "${chromosomes}"

echo $chromosomes

# Manually split the resolution string into an array
IFS=', ' read -r -a resolutions <<< "${resolutions//[\(\)]/}"

declare -A resolutionMap
for res in "${resolutions[@]}"; do
  IFS=':' read -r -a keyValue <<< "$res"
  resolutionMap[${keyValue[0]}]=${keyValue[1]}
done

# Download the GM12878 hic file if it doesn't already exist
hic_GM12878_file="../../bin/softwarefiles/GSE63525_GM12878_30.hic"
if [ ! -f "$hic_GM12878_file" ]; then
    wget $hic_GM12878_url -O "$hic_GM12878_file"
else
    echo "GM12878 hic file already exists: $hic_GM12878_file"
fi

# Download the IMR90 hic file if it doesn't already exist
hic_IMR90_file="../../bin/softwarefiles/GSE63525_IMR90_30.hic"
if [ ! -f "$hic_IMR90_file" ]; then
    wget $hic_IMR90_url -O "$hic_IMR90_file"
else
    echo "IMR90 hic file already exists: $hic_IMR90_file"
fi

path_to_jar="../../bin/softwarefiles/juicer_tools_1.22.01.jar"

# Check and run the Java commands only if the files do not already exist
for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    for chromosome in "${chromosomes[@]}"; do
        prefix="GM12878"
        path_to_hic="$hic_GM12878_file"

        eigenvector_directory="${output_directory}/eigenvector/"
        mkdir -p $eigenvector_directory
        path_to_new_eigenvector="${eigenvector_directory}res${res}_ch${chromosome}_${data_type}_${prefix}_${normalization}_eigenvector.txt"

        if [[ ! -f "$path_to_new_eigenvector" ]]; then
            eigenvector_command="java -jar ${path_to_jar} eigenvector -p ${normalization} ${path_to_hic} ${chromosome} BP ${res} ${path_to_new_eigenvector}"
            echo "Running: $eigenvector_command"
            $eigenvector_command
        else
            echo "Eigenvector file already exists: $path_to_new_eigenvector"
        fi

        prefix="IMR90"
        path_to_hic="$hic_IMR90_file"

        path_to_new_eigenvector="${eigenvector_directory}res${res}_ch${chromosome}_${data_type}_${prefix}_${normalization}_eigenvector.txt"

        if [[ ! -f "$path_to_new_eigenvector" ]]; then
            eigenvector_command="java -jar ${path_to_jar} eigenvector -p ${normalization} ${path_to_hic} ${chromosome} BP ${res} ${path_to_new_eigenvector}"
            echo "Running: $eigenvector_command"
            $eigenvector_command
        else
            echo "Eigenvector file already exists: $path_to_new_eigenvector"
        fi
    done
done

