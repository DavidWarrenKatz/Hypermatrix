#!/bin/bash
# Description: The following script pulls hi-c matrices for further processing downstream
# [TO-DO: 1]: add to readme that the getData.sh file should have chmod +x permissions 
# [TO-DO:2]: create a yml with installations for https://pypi.org/project/hic-straw/ hicstraw pip install hic-straw

# Set the path for data storage
data_path=../projects/b1198/epifluidlab/david/GSE63525/GM12878/
mkdir -p $data_path # mkdir if it doesn't exist

# Variables
resolutions=(100000)
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
data_types=("observed")
genomeID="hg19"
hic_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"

# Convert lists to strings for Python script
resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }

# Execute the Python script
python extract_hic_data.py "$data_path" "$hic_url" "$resolutions_list" "$chromosomes_list" "$data_types_list"

# Java command to create .hic files from short form text files
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        for data in "${data_types[@]}"; do
            path_to_jar="/projects/b1198/epifluidlab/david/softwareFiles/juicer_tools_1.22.01.jar"
            path_to_short_form="${data_path}hicFiles/short_score_textform/shortScore_res${resolution}_ch${chromosome}_${data}_KR.txt"
            path_to_new_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_${data}_KR.hic"
            java_command="java -Xmx100G -jar ${path_to_jar} pre -r ${resolution} ${path_to_short_form} ${path_to_new_hic} hg19"
            echo "Running: $java_command"
            $java_command
        done
    done
done

# Execute Pearson's correlation matrix
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        path_to_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_observed_KR.hic"
        path_to_new_pearsons="${data_path}pearsons/individual/res${resolution}_ch${chromosome}_observed_KR_pearsons.txt"
        pearsons_command="java -jar ${path_to_jar} pearsons -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_pearsons}"
        echo "Running: $pearsons_command"
        $pearsons_command
    done
done

# Execute eigenvector calculation
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        path_to_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_observed_KR.hic"
        path_to_new_eigenvector="${data_path}eigenvector/res${resolution}_ch${chromosome}_observed_KR_eigenvector.txt"
        eigenvector_command="java -jar ${path_to_jar} eigenvector -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_eigenvector}"
        echo "Running: $eigenvector_command"
        $eigenvector_command
    done
done

# Execute the second Python script
python process_pearsons.py "$data_path" "$resolutions_list" "$chromosomes_list" "$data_types_list"
