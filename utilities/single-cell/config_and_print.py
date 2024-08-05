# config_and_print.py
import scipy.io

bam_directory='/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
software_directory='../../bin/softwarefiles'
chrom_file="../../bin/softwarefiles/hg19.autosome.chrom.sizes"
fragments_file="$bam_directory/hg19_DpnII.txt"
output_directory='../../projects/single_cell_files'
filtered_list="$output_directory/filtered_bam_list.txt"
schicluster_env='schicluster2'
bisulfite_env='bisulfitehic27'

print(f"bam_directory='{bam_directory}'")
