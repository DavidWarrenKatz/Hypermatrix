# config_and_print.py
#This file sets the default parameters for the scNomeHic pipeline

import scipy.io

#########################################################################################################
#Define the bam directory
#The scNomeHic technique produces bam files, from which the hic and methylation information is derived. 
#The bam_directory is the ocation of these files, which must be in the form
#<prefix>.bam such as
#sc10.ACTTGA.b37.calmd.bam or
#batch1.sc15.CGATGT.hg38.calmd.bam
#The index file should be in the same directory and should be named
#<prefix>.bam.bai
#bam_directory = '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
#bam_directory = '/ocean/projects/mcb190124p/dnaase/projects/hexa_seq/scNOMe_DNA_IMR90_GM_all/bam'
#bam_directory = '/ocean/projects/mcb190124p/dnaase/startup/projects/hexa_seq/scNOMe_DNA_IMR90_GM_all/raw'
#This directory has symbolic liknks to files that I have no permisiion to read, which is why I am not using it
bam_directory = '/jet/home/dkatz/tmp_ondemand_ocean_mcb190124p_symlink/dkatz/files_from_quest/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
##########################################################################################################

##########################################################################################################
#The code needs to be updated to produce the bw methylation files from the bam files
#Right now, the pipeline assumes the methylation files are in the following directory
#and in the following form
#<prefix>.methy.bw
#such as
#batch5.scD_9.ATCACG.methy.bw
#Note, the prefixes might be a subset of the prefixes in the bam files
#methy_directory = 'methy_directory = '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam/methylation/filter_low_qual'
#The above directroy corresponds to the file which I have no access to yet
#methy_directory = '/ocean/projects/mcb190124p/dnaase/projects/hexa_seq/scNOMe_DNA_IMR90_GM_all/bam/methy'

methy_directory = '/jet/home/dkatz/tmp_ondemand_ocean_mcb190124p_symlink/dkatz/files_from_quest/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/methylation/filter_low_qual'

##########################################################################################################

software_directory = '../../bin/softwarefiles'
output_directory = '../../projects/single_cell_files'
filtered_list = f"{output_directory}/filtered_bam_list.txt"
min_high_quality_reads=250000
#resolutions = ("250000:250kb")  # Add resolutions here as a list of strings, resolution: label
resolutions = ("1000000:1Mb")
impute = False
cluster_compartments = False
cumulant = True
iterations = 400
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
mappability_threshold = 0.5
normalization = 'NONE' #options are 'NONE', 'KR', or 'ICE'
data_type = 'oe'    #options are 'oe' or 'none'
correlation = True	

#################################################################################################
#This version assumes the two cell types are GM12878 and IMR90
#This assumption needs to be updated later
hic_GM12878_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"
hic_IMR90_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5Fcombined%5F30.hic"
#################################################################################################

##########################################################################################################
#reference genome specific parameters
# Define paths based on the reference genome
reference_genome = "hg19" #needs to be either hg19 (aka GRCh37 or b37 [2009]) or hg38 (aka GRCh38 or b38 or hg20 [2013])

#The following code asusmes the DpnII restriction enzyme was used in the hic experiment. 
#For other restriction enzymes, you will need to download the restriction site file, for
#example from the juicer software files
if reference_genome == "hg19":
    chrom_file = f"{software_directory}/hg19.autosome.chrom.sizes"
    fragments_file = f"{software_directory}/hg19_DpnII.txt"
    hg_fa_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
    dark_regions_hg_url = 'https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig'
elif reference_genome == "hg38":
    chrom_file = f"{software_directory}/hg38.autosome.chrom.sizes"
    fragments_file = f"{software_directory}/hg38_DpnII.txt"
    hg_fa_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    dark_regions_hg_url = 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k36.Umap.MultiTrackMappability.bw'
else:
    raise ValueError("Invalid reference genome specified. Please choose either 'hg19' or 'hg38'.")

#########################################################################################################

#Write all parameters to a dictionary
config = {
    "bam_directory": bam_directory,
    "methy_directory": methy_directory,
    "software_directory": software_directory,
    "reference_genome": reference_genome,
    "chrom_file": chrom_file,
    "fragments_file": fragments_file,
    "output_directory": output_directory,
    "hg_fa_url": hg_fa_url,
    "filtered_list": filtered_list,
    "min_high_quality_reads": min_high_quality_reads,
    "resolutions": resolutions,
    "impute": impute,
    "cluster_compartments": cluster_compartments,
    "cumulant": cumulant,
    "iterations": iterations,
    "chromosomes": chromosomes,
    "dark_regions_hg_url": dark_regions_hg_url,
    "mappability_threshold": mappability_threshold,
    "data_type": data_type,
    "normalization": normalization,
    "correlation": correlation,
    "hic_GM12878_url": hic_GM12878_url,
    "hic_IMR90_url": hic_IMR90_url
}

#Save the dictionary for later use in MATLAB
scipy.io.savemat('config.mat', config)

#Print out the values of the dictionary
#I beleive my bash script needs the values printed out
for key, value in config.items():
    if isinstance(value, list):
        print(f"{key}={','.join(value)}")
    else:
        print(f"{key}='{value}'")


