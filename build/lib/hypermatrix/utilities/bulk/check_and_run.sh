#!/bin/bash

# Run the Python script and source the output to import the variables
eval $(python3 config_and_print.py)

# Check and run the Java commands only if the files do not already exist
for resolution in ${resolutions[@]}; do
    for chromosome in ${chromosomes[@]}; do
        for data in ${data_types[@]}; do
            hic_directory="${data_path}hicFiles/individual/"
            mkdir -p $hic_directory # mkdir if it doesn't exist
            path_to_new_hic="${hic_directory}res${resolution}_ch${chromosome}_${data}_KR.hic"

            pearson_directory="${data_path}pearsons/individual/"
            mkdir -p $pearson_directory
            path_to_new_pearsons="${pearson_directory}res${resolution}_ch${chromosome}_${data}_KR_pearsons.txt"

            eigenvector_directory="${data_path}eigenvector/"
            mkdir -p $eigenvector_directory
            path_to_new_eigenvector="${eigenvector_directory}res${resolution}_ch${chromosome}_${data}_KR_eigenvector.txt"

            if [[ ! -f "$path_to_new_hic" || ! -f "$path_to_new_pearsons" || ! -f "$path_to_new_eigenvector" ]]; then
                ./run_java_commands.sh "$data_path" "$resolution" "$chromosome" "$data"
            else
                echo "Skipping computation for Chromosome ${chromosome}, Resolution ${resolution}, Data ${data}: Files already exist."
            fi
        done
    done
done

