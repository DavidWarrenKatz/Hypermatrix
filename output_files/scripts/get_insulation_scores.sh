#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=48:00
#SBATCH --ntasks=5
#SBATCH --mem=50000
#SBATCH --error=/projects/b1198/epifluidlab/david/grant/logs/%jgetInsulationData.err
#SBATCH --output=/projects/b1198/epifluidlab/david/grant/logs/%jgetInsulationData.out

# Load the necessary modules for conda
source /etc/profile
source ~/.bashrc

# Activate the conda environment that has fanc installed
conda activate multiomics6

data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"
grant_path="/home/dwk681/workspace/grant/"
resolutions=(500000)
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
        output_file="${grant_path}Workspaces/insulation/insulation_${resolution}_ch${chromosome}_observed_KR"
        
        # Compute insulation scores
        fanc insulation "$input_file"                         "$output_file"                         -o bigwig
    done
done
