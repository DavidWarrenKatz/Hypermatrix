#!/usr/env/python 
# config_and_print.py

import argparse
import scipy.io

# Define default configuration
config = {
    "bam_directory": '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam',
    "methy_directory": '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam/methylation/filter_low_qual',
    "software_directory": '../../bin/softwarefiles',
    "chrom_file": '../../bin/softwarefiles/hg19.autosome.chrom.sizes',
    "fragments_file": '../../bin/softwarefiles/hg19_DpnII.txt',
    "output_directory": '../../projects/single_cell_files',
    "hg19_fa_url": 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
    "filtered_list": '../../projects/single_cell_files/filtered_bam_list.txt',
    "schicluster_env": 'schicluster2',
    "bisulfite_env": 'bisulfitehic27',
    "min_high_quality_reads": 250000,
    "resolutions": "1000000:1Mb",
    "impute": False,
    "cluster_compartments": False,
    "cumulant": True,
    "iterations": 400,
    "chromosomes": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'],
    "dark_regions_hg19_url": 'https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig',
    "mappability_threshold": 0.6,
    "normalization": 'NONE',
    "data_type": 'oe',
    "correlation": True,
    "genomeID": "hg19",
    "hic_GM12878_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic",
    "hic_IMR90_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5Fcombined%5F30.hic"
}

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Update config parameters via command-line arguments.')

# Add arguments for each configurable option
parser.add_argument('--methylation_file', type=str, help='Path to the single-cell CpG methylation file')
parser.add_argument('--conformation_file', type=str, help='Path to the chromosome conformation file')
parser.add_argument('--output_dir', type=str, help='Directory where the output files will be saved')
parser.add_argument('-cumulant', action='store_true', help='Run cumulant shell script')
parser.add_argument('-impute', action='store_true', help='Run impute shell script')

# Parse the arguments
args = parser.parse_args()

# Update the config with command-line arguments, if provided
if args.methylation_file:
    config["methy_directory"] = args.methylation_file
if args.conformation_file:
    config["chrom_file"] = args.conformation_file
if args.output_dir:
    config["output_directory"] = args.output_dir
if args.cumulant:
    config["cumulant"] = True
if args.impute:
    config["impute"] = True

# Save the updated config to a .mat file
scipy.io.savemat('config.mat', config)

# Print the updated config
for key, value in config.items():
    if isinstance(value, list):
        print(f"{key}={','.join(value)}")
    else:
        print(f"{key}='{value}'")



# # config_and_print.py
# import scipy.io



# bam_directory = '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
# methy_directory = '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam/methylation/filter_low_qual'
# software_directory = '../../bin/softwarefiles'
# chrom_file = f"{software_directory}/hg19.autosome.chrom.sizes"
# fragments_file = f"{software_directory}/hg19_DpnII.txt"
# output_directory = '../../projects/single_cell_files'
# hg19_fa_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
# filtered_list = f"{output_directory}/filtered_bam_list.txt"
# schicluster_env = 'schicluster2'
# bisulfite_env = 'bisulfitehic27'
# min_high_quality_reads=250000
# #resolutions = ("250000:250kb")  # Add resolutions here as a list of strings, resolution: label
# resolutions = ("1000000:1Mb")
# impute = False
# cluster_compartments = False
# cumulant = True
# iterations = 400
# chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
# dark_regions_hg19_url = 'https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig'
# mappability_threshold = 0.6
# normalization = 'NONE' #options are 'NONE', 'KR', or 'ICE'
# data_type = 'oe'    #options are 'oe' or 'none'
# correlation = True	
# genomeID = "hg19"
# hic_GM12878_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"
# hic_IMR90_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5Fcombined%5F30.hic"

# config = {
#     "bam_directory": bam_directory,
#     "methy_directory": methy_directory,
#     "software_directory": software_directory,
#     "chrom_file": chrom_file,
#     "fragments_file": fragments_file,
#     "output_directory": output_directory,
#     "hg19_fa_url": hg19_fa_url,
#     "filtered_list": filtered_list,
#     "schicluster_env": schicluster_env,
#     "bisulfite_env": bisulfite_env,
#     "min_high_quality_reads": min_high_quality_reads,
#     "resolutions": resolutions,
#     "impute": impute,
#     "cluster_compartments": cluster_compartments,
#     "cumulant": cumulant,
#     "iterations": iterations,
#     "chromosomes": chromosomes,
#     "dark_regions_hg19_url": dark_regions_hg19_url,
#     "mappability_threshold": mappability_threshold,
#     "data_type": data_type,
#     "normalization": normalization,
#     "correlation": correlation,
#     "genomeID": genomeID,
#     "hic_GM12878_url": hic_GM12878_url,
#     "hic_IMR90_url": hic_IMR90_url
# }

# scipy.io.savemat('config.mat', config)

for key, value in config.items():
    if isinstance(value, list):
        print(f"{key}={','.join(value)}")
    else:
        print(f"{key}='{value}'")


