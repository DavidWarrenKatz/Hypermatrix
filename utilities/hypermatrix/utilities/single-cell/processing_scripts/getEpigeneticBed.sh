#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=1-1:50:00
#SBATCH --ntasks=5
#SBATCH --mem=100000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetEpigeneticFromBed.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetEpigeneticFromBed.out


resolutions=(25000)
data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")
genomeID="hg19"
resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }

ftp_url="https://www.encodeproject.org/files/ENCFF682WIQ/@@download/ENCFF682WIQ.bed.gz"
hic_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119171/suppl/GSE119171_JL.merged.sorted.nodup.long.juice.q30.hic"


resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }





python << CODE
import numpy as np
import urllib.request
import hicstraw
import gzip
import io
import math
import os
import h5py
from scipy.spatial import ckdtree
from scipy.stats import ks_2samp
from multiprocessing import Pool, cpu_count
import h5py

def compute_ks_statistic(vector_i, vector_j):
    ks_statistic, p_value = ks_2samp(vector_i, vector_j)
    return ks_statistic, p_value

def calculate_ks_parallel(args):
    methylation_data, chromosome, bin_i, bin_j = args
    if not methylation_data[chromosome][f'bin{bin_i}'] or not methylation_data[chromosome][f'bin{bin_j}']:
        return bin_i, bin_j, np.nan, np.nan
    else:
        ks_statistic, p_value = compute_ks_statistic(methylation_data[chromosome][f'bin{bin_i}'], methylation_data[chromosome][f'bin{bin_j}'])
        return bin_i, bin_j, ks_statistic, p_value

def calculate_matrices_for_resolution(resolution, methylation_data, chromosomes, chrom_lengths, hic):
    path = '$data_path'
    matrices_A = {}
    matrices_B = {}
    matrices_C = {}

    for chromosome in chromosomes:
        numeric_chrom = int(chromosome[3:])
        chrom_length = chrom_lengths[chromosome]
        num_bins = math.ceil(chrom_length / resolution)
        startingbp = 0
        endingbp = chrom_length - 1
        
        matrix_filename = f"{path}Workspaces/individual/ch{numeric_chrom}_res{resolution}_oe/log2_KR.mat"
        file = h5py.File(mat_data, 'r')
        matrix_A = file['log2_matrix'][:]
        matrices_A[chromosome] = matrix_A

        num_bins = math.ceil(chrom_length / resolution)

        bin_pairs = [(i, j) for i in range(num_bins) for j in range(num_bins)]

        matrix_B = np.zeros((num_bins, num_bins))
        matrix_C = np.zeros((num_bins, num_bins))

        with Pool(cpu_count()) as pool:
            results = pool.map(calculate_ks_parallel, [(methylation_data, chromosome, bin_i, bin_j) for bin_i, bin_j in bin_pairs])

        for bin_i, bin_j, ks_statistic, p_value in results:
            matrix_B[bin_i, bin_j] = ks_statistic
            matrix_B[bin_j, bin_i] = ks_statistic
            matrix_C[bin_i, bin_j] = p_value
            matrix_C[bin_j, bin_i] = p_value

        matrices_B[chromosome] = matrix_B
        matrices_C[chromosome] = matrix_C

    return matrices_A, matrices_B, matrices_C




def main():
    path = '$data_path'
    ftp_url = '$ftp_url'
    hic_url = '$hic_url'

    # Download the methylation BED file
    local_methylation_path, _ = urllib.request.urlretrieve(ftp_url)

    # Download the Hi-C file
    hic = hicstraw.HiCFile(hic_url)
    chrom_lengths = {chrom.name: chrom.length for chrom in hic.getChromosomes()}

    resolutions = [$resolutions_list]
    chromosomes = [$chromosomes_list]
    data_types = [$data_types_list]

    for resolution in resolutions:
        resolution = int(resolution)
        matrices_A = {}
        matrices_B = {}
        matrices_C = {}

        # Read and process epigenetic data
        epigenetic_data = {}
        with gzip.open(local_methylation_path, 'rb') as gz_file:
            decompressed_data = gz_file.read().decode('utf-8')
        lines = decompressed_data.split('\n')
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                chrom, start, end, _, score, _, _, _, _, _  = parts
                if chrom in chromosomes:
                    if chrom not in epigenetic_data:
                        chrom_length = chrom_lengths[chrom]
                        num_bins = math.ceil(chrom_length / resolution)
                        epigenetic_data[chrom] = {f'bin{i}': [] for i in range(num_bins)}
                    bin_index = int(end) // resolution
                    if bin_index < len(epigenetic_data[chrom]):
                        epigenetic_data[chrom][f'bin{bin_index}'].append(float(score))

        # Calculate matrices and save
        matrices_A, matrices_B, matrices_C = calculate_matrices_for_resolution(
            resolution, epigenetic_data, chromosomes, chrom_lengths, hic)

        mdic = {"matrices_A": matrices_A, "matrices_B": matrices_B, "matrices_C": matrices_C}
        # Save matrices
        if resolution > 25000:
            savemat(os.path.join(path, f'Workspaces/pearsons_ks_pvalue_{resolution}.mat'), mdic)
            np.save(os.path.join(path, f'Workspaces/pearsons_matrices_A_{resolution}_parallel.npy'), matrices_A)
            np.save(os.path.join(path, f'Workspaces/ks_statistic_matrices_B_{resolution}_parallel.npy'), matrices_B)
            np.save(os.path.join(path, f'Workspaces/p_value_matrices_C_{resolution}_parallel.npy'), matrices_C)
        else:
            savemat(os.path.join(path, f'Workspaces/H3K9me3/H3K9me3_oe_ks_pvalue_{resolution}.mat'), mdic)         
            np.save(os.path.join(path, f'Workspaces/H3K9me3/H3K9me3_oe_matrices_A_{resolution}_parallel.npy'), matrices_A)
            np.save(os.path.join(path, f'Workspaces/H3K9me3/H3K9me3_ks_statistic_matrices_B_{resolution}_parallel.npy'), matrices_B)
            np.save(os.path.join(path, f'Workspaces/H3K9me3/H3K9me3_p_value_matrices_C_{resolution}_parallel.npy'), matrices_C)

if __name__ == "__main__":
    main()
CODE
















