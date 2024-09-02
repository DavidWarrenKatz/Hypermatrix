#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=1-01:50:00
#SBATCH --ntasks=5
#SBATCH --mem=80000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetData.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetData.out

# Set the path for data storage
data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"

resolutions=(100000)
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
data_types=("observed")
genomeID="hg19"

hic_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"

resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }

<<comment
python << EOF
import numpy as np
import hicstraw
import h5py

path = '$data_path'
hic_url = '$hic_url'
hic = hicstraw.HiCFile(hic_url)

# List of resolutions and chromosomes for which you want to create the matrices
normalization = 'KR'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
data_types = [$data_types_list]

#extract and save individual hic files
for resolution in resolutions:
    resolution = int(resolution)
    for chromosome in chromosomes:
        startingbp = 0
        endingbp = hic.getChromosomes()[int(chromosome)].length

        for data in data_types:
            matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, normalization, "BP", resolution)
            matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
            #Write short-form text file
            with open(path + f'hicFiles/short_score_textform/shortScore_res{resolution}_ch{chromosome}_{data}_{normalization}.txt', 'wt') as fid_org:
                for i in range(matrix.shape[0]):
                    for j in range(i, matrix.shape[1]):
                        if matrix[i, j] != 0:
                            fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))

            # Save the matrix as a .mat file
            new_filename = path + f'Workspaces/individual/ch{chromosome}_res{resolution}_{data}_{normalization}.mat'
            with h5py.File(new_filename, 'w') as hf:
                hf.create_dataset('matrix', data=matrix)
'''
for data in data_types:
    #extract and save a single hic file
    with open(path + f'hicFiles/short_score_textform/shortScore_all_{data}_KR.txt', 'wt') as fid_org:
        for resolution in resolutions:
            resolution = int(resolution)
            for chromosome in chromosomes:
                startingbp = 0
                endingbp = hic.getChromosomes()[int(chromosome)].length
                matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, "KR", "BP", resolution)
                matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                for i in range(matrix.shape[0]):
                    for j in range(i, matrix.shape[1]):
                        if matrix[i, j] != 0:
                            fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))
'''
EOF

#Now, run the Java command after all the short form  files are created
#for data in "${data_types[@]}"; do
#    path_to_jar="/projects/b1198/epifluidlab/david/softwareFiles/juicer_tools_1.22.01.jar"
#    path_to_short_form="${data_path}hicFiles/short_score_textform/shortScore_all_${data}_KR.txt"
#    path_to_new_hic="${data_path}hicFiles/allChromosomes_${data}.hic"
#    java_command="java -Xmx100G -jar ${path_to_jar} pre ${path_to_short_form} ${path_to_new_hic} hg19"
#    echo "Running: $java_command"
#    $java_command 
#done
comment


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


# Loop for executing the provided Java command for Pearson's correlation matrix
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        path_to_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_observed_KR.hic"
        path_to_new_pearsons="${data_path}pearsons/individual/res${resolution}_ch${chromosome}_observed_KR_pearsons.txt"
        pearsons_command="java -jar ${path_to_jar} pearsons -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_pearsons}"
        echo "Running: $pearsons_command"
        $pearsons_command
    done
done


# Loop for executing the provided Java command for eigenvector calculation
path_to_jar="/projects/b1198/epifluidlab/david/softwareFiles/juicer_tools_1.22.01.jar"
for resolution in "${resolutions[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        path_to_hic="${data_path}hicFiles/individual/res${resolution}_ch${chromosome}_observed_KR.hic"
        path_to_new_eigenvector="${data_path}eigenvector/res${resolution}_ch${chromosome}_observed_KR_eigenvector.txt"
        eigenvector_command="java -jar ${path_to_jar} eigenvector -p NONE ${path_to_hic} ${chromosome} BP ${resolution} ${path_to_new_eigenvector}"
        echo "Running: $eigenvector_command"
        $eigenvector_command
    done
done


# Start of Python code
python << PYTHON_CODE

import numpy as np
import h5py

resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
data_types = [$data_types_list]
path = '$data_path'

def read_matrix_from_text_file(file_path):
    # Read the data from the text file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Process each line and build the matrix
    matrix = []
    for line in lines:
        row = [float(val) if val != 'nan' else np.nan for val in line.strip().split()]
        matrix.append(row)

    # Convert the list of lists to a NumPy array
    matrix = np.array(matrix)

    return matrix

for chromosome in chromosomes:
    for resolution in resolutions:
        resolution = int(resolution)
        # Load the text file for each chromosome and resolution

        filename = f"{path}pearsons/individual/res{resolution}_ch{chromosome}_observed_KR_pearsons.txt"
        matrix = read_matrix_from_text_file(filename)

        translated_matrix = matrix + 1
        
        # Create a new .mat file to save the translated matrix
        new_filename = f"{path}Workspaces/individual/translated_ch{chromosome}_res{resolution}_pearsons_KR.mat"
        new_filename_npy = f"{path}Workspaces/individual/translated_ch{chromosome}_res{resolution}_pearsons_KR.npy"

        with h5py.File(new_filename, 'w') as hf:
            hf.create_dataset('translated_matrix', data=translated_matrix)

        # Save the matrix as a .npy file
        np.save(new_filename_npy, translated_matrix)
        
PYTHON_CODE




