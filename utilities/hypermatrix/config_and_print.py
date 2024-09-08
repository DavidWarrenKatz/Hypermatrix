# config_and_print.py
import scipy.io

bam_directory = '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
methy_directory = '/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam/methylation/filter_low_qual'
software_directory = '../../bin/softwarefiles'
chrom_file = f"{software_directory}/hg19.autosome.chrom.sizes"
fragments_file = f"{software_directory}/hg19_DpnII.txt"
output_directory = '../../projects/single_cell_files'
hg19_fa_url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
filtered_list = f"{output_directory}/filtered_bam_list.txt"
schicluster_env = 'schicluster2'
bisulfite_env = 'bisulfitehic27'
min_high_quality_reads=250000
#resolutions = ("250000:250kb")  # Add resolutions here as a list of strings, resolution: label
resolutions = ("1000000:1Mb")
impute = False
cluster_compartments = False
cumulant = True
iterations = 400
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
dark_regions_hg19_url = 'https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig'
mappability_threshold = 0.6
normalization = 'NONE' #options are 'NONE', 'KR', or 'ICE'
data_type = 'oe'    #options are 'oe' or 'none'
correlation = True	
genomeID = "hg19"
hic_GM12878_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"
hic_IMR90_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5Fcombined%5F30.hic"

config = {
    "bam_directory": bam_directory,
    "methy_directory": methy_directory,
    "software_directory": software_directory,
    "chrom_file": chrom_file,
    "fragments_file": fragments_file,
    "output_directory": output_directory,
    "hg19_fa_url": hg19_fa_url,
    "filtered_list": filtered_list,
    "schicluster_env": schicluster_env,
    "bisulfite_env": bisulfite_env,
    "min_high_quality_reads": min_high_quality_reads,
    "resolutions": resolutions,
    "impute": impute,
    "cluster_compartments": cluster_compartments,
    "cumulant": cumulant,
    "iterations": iterations,
    "chromosomes": chromosomes,
    "dark_regions_hg19_url": dark_regions_hg19_url,
    "mappability_threshold": mappability_threshold,
    "data_type": data_type,
    "normalization": normalization,
    "correlation": correlation,
    "genomeID": genomeID,
    "hic_GM12878_url": hic_GM12878_url,
    "hic_IMR90_url": hic_IMR90_url
}

scipy.io.savemat('config.mat', config)



for key, value in config.items():
    if isinstance(value, list):
        print(f"{key}={','.join(value)}")
    else:
        print(f"{key}='{value}'")


# To-do 
# UNCOMENT ABOVE 
# EDIT IN MAIN SO THAT IF VERBOSE ARGUMENT IS SUPPLIED THE CONFIG IS PRINTED 
# ELSE DEFAULT BEHAVIOR IS SILENT 

