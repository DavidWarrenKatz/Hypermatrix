######################################################################################                   #This shell script computes the compartment calls using the software of scHICluster
#Author: David Katz (davidkatz02@gmail.com)                                                              #####################################################################################

#!/bin/bash

eval "$(python3 config_and_print.py)"

# Source the conda environment setup script
# NEED TO UPDATE THIS SO THAT IT WORKS ON OTHER MACHINES
source /software/miniconda3/4.12.0/etc/profile.d/conda.sh

schicluster_env=schicluster2
bisulfite_env=bisulfitehic27

conda activate $bisulfite_env

# Capture the Python path
bisulfite_env_python_path=$(which python)
echo "Python path in bisulfite_env: $bisulfite_env_python_path"

#########################################################################################                ###compartment calling                                                                                   #########################################################################################

# First, check if the hg19.fa.gz file exists in the software directory
if [ ! -f "$software_directory/hg19.fa.gz" ]; then
  echo "Downloading hg19.fa.gz into the software directory..."
  wget --no-passive-ftp -P $software_directory $hg19_fa_url 
else
  echo "hg19.fa.gz already exists in the software directory."
fi

# Uncompress the hg19.fa.gz file if it is not already uncompressed
if [ ! -f "$software_directory/hg19.fa" ]; then
  echo "Uncompressing hg19.fa.gz..."
  gunzip "$software_directory/hg19.fa.gz"
else
  echo "hg19.fa already exists in the software directory."
fi

# Path to the downloaded and uncompressed hg19 reference genome
hg19_fa_path="$software_directory/hg19.fa"

# I beleive the scHiCluster code recommends using the imputed matrix, but I am not doing that here
# Split the resolutions string into an array and iterate over it
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    # Split each pair into resolution and label
    IFS=':' read -r resolution label <<< "$pair"
    # Create the output directory variable for the current resolution
    if [ "$impute" = "True" ]; then
        dir="./hicluster_${label}_impute_dir/bins"
    else
        dir="./hic_${label}_raw_dir/bins"
    fi
    
    # Create the directory
    mkdir -p "$dir"
    echo "Created directory: $dir"    
done

# Execute the Python script                                                                                            
python process_chromsizes.py $dir

# Split the resolutions string into an array and iterate over it
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    # Split each pair into resolution and label
    IFS=':' read -r resolution label <<< "$pair"
    # Concatenate and sort all chromosome BED files into a single file for the current resolution
    # Sorting is based on chromosome number and start position
    cat ./hicluster_${label}_impute_dir/bins/chr*.bedi | sort -k1,1 -k2,2n > ./hicluster_${label}_impute_dir/bins/hg19.${label}_bin.bed
    
    # Use bedtools to calculate CG density for the segments defined in the sorted BED file
    # This requires a reference genome file and the sorted BED file as inputs
    bedtools nuc -fi $software_directory/genomes/ucsc_hg19/hg19.fa -bed $output_directory/hicluster_${label}_impute_dir/bins/hg19.${label}_bin.bed -pattern CG -C > $output_directory/hicluster_${label}_impute_dir/bins/hg19.${label}_bin.cg_density.bed
done

# Execute the hicluster comp-cpg-cell command
ls sc*.hic_matrix.txt.gz | perl -sne '
  use strict;
  use warnings;

  # Accept resolutions from command line arguments
  our $resolutions_string;
  our $software_directory;
  our $submission_script_directory;
  our $schicluster_env;
  
  # Convert resolutions string to a hash
  my %resolutions = map { split /:/ } split /,/, $resolutions_string;

  # Remove newline and store the input filename
  chomp(my $in = $_);

  # Extract base ID by removing the file extension
  (my $id = $in) =~ s/\.hic_matrix\.txt\.gz$//;

  # Loop over chromosomes and resolutions
  my @chromosomes = (1..22);
  foreach my $res (keys %resolutions) {
      my $label = $resolutions{$res};
      my $indir = "hicluster_${label}_impute_dir/";
      my $outdir = $indir;
      my $cpg_file = "${indir}bins/hg19.${label}_bin.cg_density.bed";
      foreach my $c (@chromosomes) {
          # Construct a sanitized job name to avoid issues with special characters
          my $job_name = "${label}_impute_${id}.${c}";
          $job_name =~ s/[^\w]/_/g;

          # Build the command for Hi-C CPG density computation
          my $cmd = "source activate $schicluster_env && hicluster comp-cpg-cell ".
                    "--indir $indir ".
                    "--outdir $outdir ".
                    "--cell $id ".
                    "--chrom chr$c ".
                    "--mode pad1_std1_rp0.5_sqrtvc ".
                    "--cpg_file $cpg_file";

          # Submit the job through the ssubexc.pl script with sanitized job name
          my $submit_cmd = "perl $submission_script_directory/ssubexc.pl ".
                           "\"$cmd\" \"$job_name\" 4 2";
          system($submit_cmd);
      }
  }
' -- -resolutions_string="$resolutions" -software_directory="$software_directory" -submission_script_directory="$submission_script_directory" -schicluster_env="$schicluster_env"

# Split the resolutions string into an array and iterate over it
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    # Split each pair into resolution and label
    IFS=':' read -r resolution label <<< "$pair"
    
    # Create the filelist and merged directories for each resolution
    mkdir -p "hicluster_${label}_impute_dir/filelist"
    mkdir -p "hicluster_${label}_impute_dir/merged"

    conda activate schicluster2

    # Using Perl to generate file lists for each chromosome for the current resolution
    perl -e "\@cs=(1..22);foreach \$c(\@cs){\`ls hicluster_${label}_impute_dir/chr\${c}/*_chr\${c}_pad1_std1_rp0.5_sqrtvc.cpgcomp.npy > hicluster_${label}_impute_dir/filelist/complist_pad1_std1_rp0.5_sqrtvc_chr\$c.txt\`}"

    # Concatenate cell lists for each chromosome for the current resolution
    perl -e "\@cs=(1..22);foreach \$c(\@cs){print \"\$c\n\";\`hicluster comp-concatcell-chr --cell_list hicluster_${label}_impute_dir/filelist/complist_pad1_std1_rp0.5_sqrtvc_chr\$c.txt --outprefix hicluster_${label}_impute_dir/merged/pad1_std1_rp0.5_sqrtvc_chr\$c --ncpus 10\`;}"

done

# Create the Python script that will make a merged file
cat > make_merged.py <<EOF
import numpy as np
import os
import argparse
from sklearn.preprocessing import quantile_transform

# Parse command-line arguments for resolution and label
parser = argparse.ArgumentParser(description='Merge compartment calls for different resolutions.')
parser.add_argument('--resolution', type=str, help='Resolution value, e.g., "1000000" for 1Mb')
parser.add_argument('--label', type=str, help='Resolution label, e.g., "1Mb"')
parser.add_argument('--indir_base', type=str, help='Base input directory, e.g., "/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/new_processing3/"')

args = parser.parse_args()

# Set variables based on command-line arguments
res0 = args.label
mode = 'pad1_std1_rp0.5_sqrtvc'
indir = f'{args.indir_base}hicluster_{res0}_impute_dir/merged/'

comp = []
for c in range(1, 23):
    tmp_file = f'{indir}{mode}_chr{c}.cpgcomp.npy'
    if os.path.exists(tmp_file):
        tmp = np.load(tmp_file)
        binfilter = (np.std(tmp, axis=0) > 0)
        comptmp = np.ones(tmp.shape) / 2
        comptmp[:, binfilter] = quantile_transform(tmp[:, binfilter], output_distribution='uniform', n_quantiles=int(np.sum(binfilter)//10), axis=1)
        comp.append(comptmp)
        print(c)

comp = np.concatenate(comp, axis=1)
output_file = f'{indir}all_merged_{res0}_cpgcmp.{mode}.txt'
np.savetxt(output_file, comp, fmt='%s', delimiter='\t')
EOF

# Iterate over the resolutions to execute the Python script
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    IFS=':' read -r resolution label <<< "$pair"
    
    # Execute the Python script for the current resolution
    python make_merged.py --resolution "${resolution}" --label "${label}" --indir_base "${indir_base}"
done


# Create the Python script that will save the merged matrix
cat > save_merged_matrix.py <<EOF
import os
import glob
import numpy as np
import argparse
import re

# Parse command-line arguments for resolution, label, and base directory
parser = argparse.ArgumentParser(description='Save merged matrix for different resolutions.')
parser.add_argument('--label', type=str, help='Resolution label, e.g., "1Mb"')
parser.add_argument('--base_dir', type=str, help='Base directory for input and output files')
args = parser.parse_args()

label = args.label
base_dir = os.path.join(args.base_dir, f'hicluster_{label}_impute_dir/merged/')

# Pattern to match all .npy files in the merged directory for the specified resolution
pattern = os.path.join(base_dir, f'pad1_std1_rp0.5_sqrtvc_chr*.cpgcomp.npy')

# Find all .npy files matching the pattern
npy_files = glob.glob(pattern)

# Function to extract the chromosome number for sorting
def extract_chr_num(filename):
    match = re.search(r'chr(\d+|\X)', filename)
    if match:
        chr_num = match.group(1)
        return int(chr_num) if chr_num.isdigit() else 25 if chr_num == 'X' else 0
    return 0  # Default case for sorting

# Sort files based on the chromosome number
npy_files_sorted = sorted(npy_files, key=extract_chr_num)

# Load each .npy file, drop the last bin, and append to the list
loaded_arrays = [np.delete(np.load(file), -1, axis=1) for file in npy_files_sorted]

# Concatenate all arrays along the first axis (rows) to match desired dimensions
large_matrix = np.concatenate(loaded_arrays, axis=1)
large_matrix = large_matrix.T

# Verify the shape of the large_matrix
print(f"The shape of the large matrix is: {large_matrix.shape}")

# Specify the output file path for saving the matrix, adjusted for the specified resolution
output_dir = os.path.join(args.base_dir, 'matrices', label)
os.makedirs(output_dir, exist_ok=True)
output_file_path = os.path.join(output_dir, f'b37.autosome.{label}_interval.add_value.hic.GM_IMR90_153_samples.npy')

# Save the matrix to a file
np.save(output_file_path, large_matrix)
EOF

# Iterate over the resolutions to execute the Python script
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    IFS=':' read -r resolution label <<< "$pair"
    
    # Execute the Python script for the current resolution
    python save_merged_matrix.py --label "${label}" --base_dir "${base_dir}"
done


