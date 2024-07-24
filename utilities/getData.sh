#!/bin/bash
# Description: The following script pulls hi-c matrices for further processing downstream

# Activate conda environment 
# eval $"(conda shell.bash hook)"
# conda activate hypermatrix

# Variables
# Run the Python script and source the output to import the variables
eval $(python3 config_and_print.py)

# Set the path for data storage
mkdir -p $data_path # mkdir if it doesn't exist

# Convert lists to strings for Python script
resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }

short_form_directory="${data_path}hicFiles/short_score_textform/"
mkdir -p $short_form_directory # mkdir if it doesn't exist
            
individual_workspace_directory="${data_path}Workspaces/individual/"
mkdir -p $individual_workspace_directory # mkdir if it doesn't exist


# Download the dark regions file if it doesn't already exist
dark_regions_file="../projects/softwarefiles/ENCFF000EHJ_hg19_wgEncodeCrgMapabilityAlign36mer.bigWig"
if [ ! -f "$dark_regions_file" ]; then
    wget https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig -O "$dark_regions_file"
fi

# Execute the Python script to create dark regions file
python retrieve_dark_bins.py

# Execute the Python script
python extract_hic_data.py "$data_path" "$hic_url" "$resolutions_list" "$chromosomes_list" "$data_types_list"

# Java command to create .hic files from short form text files
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        for data in "${data_types[@]}"; do
            path_to_jar="../projects/softwarefiles/juicer_tools_1.22.01.jar"
            short_form_directory="${data_path}hicFiles/short_score_textform/"
            path_to_short_form="${short_form_directory}shortScore_res${resolution}_ch${chromosome}_${data}_KR.txt"
            hic_directory="${data_path}hicFiles/individual/"
            mkdir -p $hic_directory # mkdir if it doesn't exist
            path_to_new_hic="${hic_directory}res${resolution}_ch${chromosome}_${data}_KR.hic"
            java_command="java -Xmx100G -jar ${path_to_jar} pre -r ${resolution} ${path_to_short_form} ${path_to_new_hic} hg19"
            echo "Running: $java_command"
            $java_command
        done
    done
done

# Execute Pearson's correlation matrix
echo -e "[LOG]: Running Pearson's correlation matrix "
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        for data in "${data_types[@]}"; do
            path_to_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_${data}_KR.hic"
            pearson_directory="${data_path}pearsons/individual/"
            mkdir -p $pearson_directory
            path_to_new_pearsons="${pearson_directory}res${resolution}_ch${chromosome}_${data}_KR_pearsons.txt"
            pearsons_command="java -jar ${path_to_jar} pearsons -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_pearsons}"
            echo "Running: $pearsons_command"
            $pearsons_command
        done
    done
done

# Execute eigenvector calculation
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        for data in "${data_types[@]}"; do
            path_to_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_${data}_KR.hic"
            eigenvector_directory="${data_path}eigenvector/"
            mkdir -p $eigenvector_directory
            path_to_new_eigenvector="${eigenvector_directory}res${resolution}_ch${chromosome}_${data}_KR_eigenvector.txt"
            eigenvector_command="java -jar ${path_to_jar} eigenvector -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_eigenvector}"
            echo "Running: $eigenvector_command"
            $eigenvector_command
        done
    done
done

# Execute the second Python script
python process_pearson.py "$data_path" "$resolutions_list" "$chromosomes_list" "$data_types_list"
