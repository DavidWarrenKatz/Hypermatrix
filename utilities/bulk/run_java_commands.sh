#!/bin/bash

data_path=$1
resolution=$2
chromosome=$3
data=$4

path_to_jar="../../bin/softwarefiles/juicer_tools_1.22.01.jar"
short_form_directory="${data_path}hicFiles/short_score_textform/"
path_to_short_form="${short_form_directory}shortScore_res${resolution}_ch${chromosome}_${data}_KR.txt"
hic_directory="${data_path}hicFiles/individual/"
path_to_new_hic="${hic_directory}res${resolution}_ch${chromosome}_${data}_KR.hic"

# Create .hic files from short form text files
if [[ ! -f "$path_to_new_hic" ]]; then
    java_command="java -Xmx100G -jar ${path_to_jar} pre -r ${resolution} ${path_to_short_form} ${path_to_new_hic} hg19"
    echo "Running: $java_command"
    $java_command
else
    echo "HIC file already exists: $path_to_new_hic"
fi

# Execute Pearson's correlation matrix
pearson_directory="${data_path}pearsons/individual/"
path_to_new_pearsons="${pearson_directory}res${resolution}_ch${chromosome}_${data}_KR_pearsons.txt"

if [[ ! -f "$path_to_new_pearsons" ]]; then
    pearsons_command="java -jar ${path_to_jar} pearsons -p NONE ${path_to_new_hic} ${chromosome} BP ${resolution} ${path_to_new_pearsons}"
    echo "Running: $pearsons_command"
    $pearsons_command
else
    echo "Pearson's file already exists: $path_to_new_pearsons"
fi

# Execute eigenvector calculation
eigenvector_directory="${data_path}eigenvector/"
path_to_new_eigenvector="${eigenvector_directory}res${resolution}_ch${chromosome}_${data}_KR_eigenvector.txt"

if [[ ! -f "$path_to_new_eigenvector" ]]; then
    eigenvector_command="java -jar ${path_to_jar} eigenvector -p NONE ${path_to_new_hic} ${chromosome} BP ${resolution} ${path_to_new_eigenvector}"
    echo "Running: $eigenvector_command"
    $eigenvector_command
else
    echo "Eigenvector file already exists: $path_to_new_eigenvector"
fi

