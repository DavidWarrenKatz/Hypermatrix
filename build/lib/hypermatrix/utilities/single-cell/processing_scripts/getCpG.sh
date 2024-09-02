#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=1-1:50:00
#SBATCH --ntasks=5
#SBATCH --mem=100000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetPearsons.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetPearsons.out


data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
genomeID="hg19"
resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }

ftp_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119171/suppl/GSE119171%5FF123%2EMethylHic%2Emm9%2Ecalmd%2Ecpg%2Efiltered%2Esort%2ECG%2Estrand%2E6plus2%2Ebed%2Egz"
hic_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119171/suppl/GSE119171_JL.merged.sorted.nodup.long.juice.q30.hic"

resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }

# Python code to import Hic_Data both HighRes and LowRes
python << EOF
import numpy as np
import hicstraw
from scipy.io import savemat

path = '$data_path'
ftp_url = '$ftp_url'
hic_url = '$hic_url'

hic = hicstraw.HiCFile(hic_url)

# List of resolutions and chromosomes for which you want to create the matrices
normalization = 'NONE'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
data_types = [$data_types_list]


for resolution in resolutions:
    resolution = int(resolution)
    for chromosome in chromosomes:
        numeric_chrom = int(chromosome[3:])
        startingbp = 0
        endingbp = hic.getChromosomes()[numeric_chrom].length

        for data in data_types:
            if resolution < 25001:
                matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, "NONE", "BP", resolution)
                large_matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                # Calculate the number of smaller blocks needed
                block_size = 500
                num_blocks_x = int(np.ceil(large_matrix.shape[0] / block_size))
                num_blocks_y = int(np.ceil(large_matrix.shape[1] / block_size))
                #Loop through the smaller blocks, saving each one
                for i in range(num_blocks_x):
                    for j in range(num_blocks_y):
                        start_i = i * block_size
                        end_i = min(start_i + block_size, large_matrix.shape[0])
                        start_j = j * block_size
                        end_j = min(start_j + block_size, large_matrix.shape[1])
                        block = large_matrix[start_i:end_i, start_j:end_j]
                        savemat(path + f'Workspaces/submatrices/ch{numeric_chrom}_res{resolution}_{data}_NONE_block_{i}_{j}.mat', {'block': block})
                savemat(path + f'Workspaces/submatrices/ch{numeric_chrom}_res{resolution}_{data}_NONE_large_matrix_size.mat', {'large_matrix_size': large_matrix.shape[0], 'block_size': block_size}  )         
            else:
                matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, "NONE", "BP", resolution)
                matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                # Write short-form text file
                with open(f'hicFiles/short_score_textform/shortScore_res{resolution}_ch{numeric_chrom}_{data}.txt', 'wt') as fid_org:
                    for i in range(matrix.shape[0]):
                        for j in range(i, matrix.shape[1]):
                            if matrix[i, j] != 0:
                                fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))
                # Save the matrix as a .mat file
                mdic = {'matrix': matrix}
                                                                                                                                                1,1           Top

