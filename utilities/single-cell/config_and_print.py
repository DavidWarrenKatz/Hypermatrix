# config_and_print.py
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
resolutions = ("1000000:1Mb")  # Add resolutions here as a list of strings, resolution: label
impute = True
cluster_compartments = False
cumulant = False

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
}


for key, value in config.items():
    if isinstance(value, list):
        print(f"{key}={','.join(value)}")
    else:
        print(f"{key}='{value}'")


