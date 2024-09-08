#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=48:00
#SBATCH --ntasks=5
#SBATCH --mem=50000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%jgetInsulationData.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%jgetInsulationData.out

data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"

resolutions=(1000000)
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
data_types=("observed")
genomeID="hg19"

resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }

for chromosome in "${chromosomes[@]}"; do
    for resolution in "${resolutions[@]}"; do
        filename="hicFiles/individual/res${resolution}_ch${chromosome}_observed_KR"
        input_file="${data_path}${filename}.hic"
        output_file="${data_path}insulation/insulation_${resolution}_ch${chromosome}_observed_KR"
        
        # Compute insulation scores
        fanc insulation "$input_file" \
                        "$output_file" \
                        -o bigwig
    done
done


<<comment
python << FILE
import pyBigWig
import numpy as np
from scipy.io import savemat
import math

data_path = '$data_path'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
#insulation_windows = ["15mb", "10mb", "7mb", "5mb", "3mb"]  #for 25,000
insulation_windows = ["1.5mb", "2.5mb", "3.5mb", "5mb", "7.5mb"] #for500,000

bin_averages_by_ch = {chromosome: {insulation_window: [] for insulation_window in insulation_windows} for chromosome in chromosomes}

def calculate_bin_averages(data, elements_per_bin):
    num_bins = math.ceil(len(data)/elements_per_bin)
    bin_averages = np.zeros(num_bins)
    for i in range(num_bins):
        start_index = i * elements_per_bin
        end_index = min((i + 1) * elements_per_bin, len(data))
        bin_data = data[start_index:end_index]
        if len(bin_data) > 0:
            bin_averages[i] = np.mean(bin_data)
        else:
            bin_averages[i] = 0
    return np.nan_to_num(bin_averages, nan=0.0)

for resolution in resolutions:
    resolution = int(resolution)
    for chromosome in chromosomes:
        mdic = {}
        for insulation_window in insulation_windows:
            bigwigfile = data_path + f'insulation/insulation_{resolution}_ch{chromosome}_observed_KR_{insulation_window}.bigwig'
            bw_file = pyBigWig.open(bigwigfile)
            chrom_sizes = bw_file.chroms()
            
            if chromosome_name in chrom_sizes:
                start_position = 0
                end_position = chrom_sizes[chromosome]
                values = np.array(bw_file.values(chromosome, start_position, end_position))
                bin_averages_by_ch[chromosome][insulation_window] = calculate_bin_averages(values, resolution) 

    # Gather all variables you want to save into a single dictionary
    data_to_save = {
    'resolution': resolution,
    'bin_averages_data_by_ch': bin_averages_data_by_ch
    }
    # Save the data to an npy file
    np.save(path + f'Images/InsulationData_res{resolution}_observed_KR.npy', data_to_save)     
FILE
comment








