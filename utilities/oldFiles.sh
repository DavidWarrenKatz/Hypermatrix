va command to create .hic files from short form text files
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
