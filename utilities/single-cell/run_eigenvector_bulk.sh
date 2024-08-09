#!/bin/bash

data_path=$1
resolution=$2
chromosome=$3
data=$4
prefix=$5

path_to_jar="../projects/softwarefiles/juicer_tools_1.22.01.jar"
path_to_new_hic="../../bin/softwarefiles/GSE63525_${prefix}IMR90.hic"
eigenvector_directory="${output_directory}eigenvector/"
path_to_new_eigenvector="${eigenvector_directory}res${resolution}_ch${chromosome}_${data}_${prefix}_KR_eigenvector.txt"

if [[ ! -f "$path_to_new_eigenvector" ]]; then
    eigenvector_command="java -jar ${path_to_jar} eigenvector -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_eigenvector}"
    echo "Running: $eigenvector_command"
    $eigenvector_command
else
    echo "Eigenvector file already exists: $path_to_new_eigenvector"
fi
