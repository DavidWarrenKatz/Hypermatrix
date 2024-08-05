######################################################################################
#This shell script performs the scHiCluster pipeline in single cell files
#For various resolutions
#Author: David Katz (davidkatz02@gmail.com)
#modified from Yaping Liu's code 
#####################################################################################

#This script assumes the hic files are in the form sc10.ACTTGA.b37.calmd.bam
#scD6.TAGCTT.b37.calmd.bam 
#This script assumes I have the restriction file hg19_DpnII.txt in the bam directory
#This script assumes I have the conda environments bisulfitehic and schicluster

#########################################################################################
###make submission script and define global variables
#This section contains all the machine-dependnet variables, that is the variables that you need to edit
########################################################################################

schicluster_env=schicluster2
bisulfite_env=bisulfitehic27

conda activate $bisulfite_env

# Capture the Python path
bisulfite_env_python_path=$(which python)
echo "Python path in bisulfite_env: $bisulfite_env_python_path"

bam_directory='/projects/b1198/epifluidlab/david/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/processed_new3'
filtered_list='/projects/b1198/epifluidlab/david/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/processed_new3/name_order_and_coloring.txt'
software_directory='/home/dwk681/workspace/softwareFiles'
submission_script_directory='.'
chrom_file="/home/dwk681/workspace/Genomic_Files/ucsc_hg19/hg19.autosome.chrom.sizes"

# Define the path where the Perl script will be created
perl_script_name="ssubexc.pl"

# Create the Perl script using a heredoc
cat <<'EOF' > "$perl_script_name"
my $cmd=$ARGV[0];
my $name=$ARGV[1];

my $random=int(rand(100000000));
my $error_log="/projects/b1198/epifluidlab/david/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/processed_new/logs/${random}ssubexc.err";
my $output_log="/projects/b1198/epifluidlab/david/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/processed_new/logs/${random}ssubexc.out";

my $PWD=`pwd`;
chomp($PWD);

my $header = "#!/usr/bin/env bash\n";
$header .= "#SBATCH --account=b1042\n";
$header .= "#SBATCH --partition=genomics\n";
$header .= "#SBATCH --time=2:00:00\n";
$header .= "#SBATCH --ntasks=5\n";
$header .= "#SBATCH --mem=19000\n";
$header .= "#SBATCH --error=$error_log\n";
$header .= "#SBATCH --output=$output_log\n";
$header .= "#SBATCH -J $name\n";
$header .= "#SBATCH --export=ALL\n";
$header .= "#SBATCH --requeue\n";

$header .= "source ~/.bash_profile\n";
$header .= "$cmd\n";

print STDERR "$cmd\n";

my $sub_file="ssub_".$name.".$random.sh";
open(O,">$sub_file") or die "can not open file:$!";
print O "$header\n";
close(O);

system("sbatch $sub_file")==0 || die "Job submission failed\n";

print STDERR "job is submitted\n";
system("unlink $sub_file");
EOF

# Make the Perl script executable
chmod +x "$perl_script_name"

echo "Perl script '$perl_script_name' has been created and made executable."

# Create preprocess.sh with the specified content
cat << 'EOF' > preprocess.sh
#!/bin/bash

# Preprocessing command
cut -d' ' -f1-7 $1 | sed 's/ /\t/g' | perl -ne 'chomp;@f=split "\t";$f[1]="chr$f[1]";$f[5]="chr$f[5]";print join("\t",@f)."\n";' | gzip -c > $2
EOF

# Make preprocess.sh executable
chmod +x preprocess.sh

# Debug: Print the path to ensure it's captured correctly
echo "Python path in bisulfite_env: $bisulfite_env_python_path"
echo "preprocess.sh has been created and made executable."

#########################################################################################
###generate contact matrix file for each single-cell from bam file
###this block is reolusiton independent, resolutions are specified in the next block
########################################################################################
mkdir logs

#create symbolic link to bam files in current directory
#include only the bam files that are listed in filtered list
# list all the files in the bam directroy and execute a Perl script creating symbolic links
ls $bam_directory/sc*.b37.calmd.bam | perl -ne '
    BEGIN {
        # Get the filtered list file path from the command line arguments
        $filtered_list = shift @ARGV;
        # Open the filtered list file or die with an error message
        open(my $fh, "<", $filtered_list) or die "Could not open file '\''$filtered_list'\'' $!";
        # Read the filtered list into a hash for quick lookups
        %dict = map { chomp; my @f = split(/\s+/); $f[1] => 1 } <$fh>;
        close($fh);
    }
    # Process each line of input from ls
    chomp;
    my $in = $_;
    # Extract the identifier from the file name
    if ($in =~ /.*\/(sc\d+\.\w+)\.b37\.calmd\.bam$/) {
        my $p = $1;
        # If the identifier is in the filtered list, create a symbolic link
        if (exists $dict{$p}) {
            my $cmd = "ln -s $in ./$p.bam";
            system($cmd);
        }
    }
' $filtered_list

# Iterate over BAM files matching the pattern "sc*.bam"
# submit a job to convert each file to the juicer_short format
ls sc*.bam | perl -sne '
  # Get variables from command line arguments
  my $software_directory = $software_directory;
  my $bam_directory = $bam_directory;
  my $filtered_list = $filtered_list;
  my $submission_script_directory = $submission_script_directory;
  my $schicluster_env = $schicluster_env;
  my $bisulfite_env_python_path = $bisulfite_env_python_path;

  chomp(my $in = $_); # Remove newline, store the original filename

  # Load samtools module before executing commands that require samtools
  my $load_samtools = "module load samtools";

  # Extract prefix from filename for job naming and command construction
  (my $p = $in) =~ s/\.bam$//;

  # Construct the command to be executed, including loading samtools module
  my $cmd = "$load_samtools && samtools view --threads 5 -bh -q 30 -f 1 -F 1804 $in > ${p}.good_reads.bam && $bisulfite_env_python_path $software_directory/sam2juicer.py -s ${p}.good_reads.bam -f $bam_directory/hg19_DpnII.txt > ${p}.hic.txt";

  # Job name construction to avoid syntax errors and ensure uniqueness
  my $job_name = "sam2juicer_short_$p";

  # Submit the job through the ssubexc.pl script
  my $submit_cmd = "perl $submission_script_directory/ssubexc.pl \"$cmd\" \"$job_name\"";
  system($submit_cmd);
' -- -software_directory=$software_directory -bam_directory=$bam_directory -filtered_list=$filtered_list -submission_script_directory=$submission_script_directory -schicluster_env=$schicluster_env -bisulfite_env_python_path=$bisulfite_env_python_path

echo "Waiting for 40 minutes before running the checks..."
sleep 2400 # Wait for 2400 seconds (40 minutes)

# Initialize counters
successful_count=0
unsuccessful_files=()

# Iterate over the base names from the filtered_list
while IFS=$'\t' read -r idx base_name color; do
  # Append .bam extension to check for the files
  bam_file="$base_name.bam"
  good_reads_file="$base_name.good_reads.bam"
  hic_txt_file="$base_name.hic.txt"

  # Check if the files exist and are not empty
  if [[ -f "$bam_file" ]] && [[ -s "$bam_file" ]] && [[ -f "$good_reads_file" ]] && [[ -s "$good_reads_file" ]] && [[ -f "$hic_txt_file" ]] && [[ -s "$hic_txt_file" ]]; then
    ((successful_count++))
  else
    unsuccessful_files+=("$base_name")
  fi
done < "$filtered_list"

# Print the results
echo "Number of successfully processed bam files: $successful_count"
if [ ${#unsuccessful_files[@]} -eq 0 ]; then
  echo "All files were processed successfully."
else
  echo "Files not processed successfully:"
  for file in "${unsuccessful_files[@]}"; do
    echo "$file"
  done
fi

# Iterate over files matching the pattern "sc*.hic.txt"
#The jobs submitted in the previous step need to be finished before this step is run
#final preprocesses to get into a format for scHiCluster
ls sc*.hic.txt | perl -sne '
  # Get variables from command line arguments
  my $software_directory = $software_directory;
  my $submission_script_directory = $submission_script_directory;

  chomp(my $in = $_); # Remove newline and store the input filename

  # Prepare the output filename by replacing ".hic.txt" with ".hic_matrix.txt.gz"
  (my $o = $in) =~ s/\.hic.txt$/.hic_matrix.txt.gz/;

  # Construct a sanitized job name to avoid issues with special characters
  (my $job_name = $in) =~ s/[^\w]/_/g;
  $job_name = "process_" . $job_name;

  # Command to preprocess the .hic.txt file into a Hi-C matrix and compress it
  my $cmd = "./preprocess.sh $in $o";

  # Submit the job through the ssubexc.pl script with sanitized job name
  my $submit_cmd = "perl $submission_script_directory/ssubexc.pl \"$cmd\" \"$job_name\" 2 1";
  system($submit_cmd);
' -- -software_directory=$software_directory -submission_script_directory=$submission_script_directory

#########################################################################################
###hicluster generatematrix-cell
###adjusts the resolution parameters as needed here, following the default format
#########################################################################################

#This is a string of resolutions and names fot the directory of each resolution
resolutions="1000000:1Mb,100000:100kb,10000:10kb"
resolutions="250000:250kb"

echo "$resolutions"
#The command generatematrix-cell makes the directroy hicluster_1Mb_raw_dir and a subdirectory for each chromosome. 
#It then populates these directories with the contact matrices in txt format. A single command produces the 
#contact matrix for each chromosome. 

ls sc*.hic_matrix.txt.gz | perl -sne '
  # Parse the resolutions string into a hash
  my %resolutions = map { split /:/ } split /,/, $resolutions_string;

  # Other script contents remain the same...
  # Remove newline and store the input filename
  chomp(my $in = $_);

  # Extract base ID by removing the file extension
  (my $id = $in) =~ s/\.hic_matrix\.txt\.gz$//;

  foreach my $res (keys %resolutions) {
      my $label = $resolutions{$res};

      # Construct a sanitized job name and other logic follows...
      (my $job_name = $id) =~ s/[^\w]/_/g;
      $job_name = "${label}_gene_mat_" . $job_name;

      my $cmd = "source activate $schicluster_env && hicluster generatematrix-cell ".
                "--infile $in ".
                "--outdir hicluster_${label}_raw_dir/ ".
                "--chrom_file $chrom_file ".
                "--res $res ".
                "--cell $id ".
                "--chr1 1 ".
                "--pos1 2 ".
                "--chr2 5 ".
                "--pos2 6";

      my $submit_cmd = "perl $submission_script_directory/ssubexc.pl ".
                       "\"$cmd\" \"$job_name\" 10 5";
      system($submit_cmd);
  }
' -- -resolutions_string="$resolutions" -software_directory="$software_directory" -submission_script_directory="$submission_script_directory" -chrom_file="$chrom_file" -schicluster_env="$schicluster_env"

# Convert resolutions string to an associative array
declare -A resolutionMap
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    IFS=':' read -r -a keyValue <<< "$pair"
    resolutionMap[${keyValue[0]}]=${keyValue[1]}
done

# Initialize counters and prepare to track unsuccessful basenames
declare -A successful_count
declare -A unsuccessful_basenames

for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    outdir="hicluster_${label}_raw_dir"
    successful_count[$res]=0
    unsuccessful_key="${res}_unsuccessful"

    # Read each line from the filtered list
    while IFS=$'\t' read -r idx base_name color; do
        expected_file="$outdir/chr1/${base_name}_chr1.txt"

        # Check if the expected file exists and is not empty
        if [[ -f "$expected_file" ]] && [[ -s "$expected_file" ]]; then
            ((successful_count[$res]++))
        else
            # Append unsuccessful basename to string associated with this resolution
            unsuccessful_basenames[$unsuccessful_key]+="$base_name "
            echo "$base_name is unsuccessful at resolution $label"
            # Optional: Attempt to reprocess unsuccessful basenames here
            # Build the command to reprocess the file
            cmd="hicluster generatematrix-cell --infile ${base_name}.hic_matrix.txt.gz --outdir $outdir/ --chrom_file $chrom_file --res $res --cell ${base_name}.hic.txt --chr1 1 --pos1 2 --chr2 5 --pos2 6"

            # Submit the job for reprocessing
            submit_cmd="perl $submission_script_directory/ssubexc.pl \"$cmd\" \"process_${base_name}_${label}\" 10 5"
            # Uncomment the line below to actually execute the command
            eval $submit_cmd
        fi
    done < "$filtered_list"

    # Print the results for the current resolution
    echo "Number of successfully processed basenames for $label: ${successful_count[$res]}"
    if [ -z "${unsuccessful_basenames[$unsuccessful_key]}" ]; then
        echo "All basenames were processed successfully at resolution $label."
    else
        echo "Basenames not processed successfully at resolution $label:"
        # Split the string back into an array for listing
        IFS=' ' read -r -a unsuccessful_array <<< "${unsuccessful_basenames[$unsuccessful_key]}"
        for base_name in "${unsuccessful_array[@]}"; do
            echo "$base_name"
        done
    fi
done

echo "Waiting for 40 minutes before running the checks..."
# Adjusted to 40 minutes as per the message
sleep 2400 # Wait for 2400 seconds (40 minutes)

# Convert resolutions string to an associative array
declare -A resolutionMap
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    IFS=':' read -r -a keyValue <<< "$pair"
    resolutionMap[${keyValue[0]}]=${keyValue[1]}
done

# List of chromosomes
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

# Iterate through resolutions
for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    outdir="hicluster_${label}_raw_dir" # Directory name includes resolution label

    # Initialize counters and arrays for this resolution
    declare -A file_creation_count
    declare -A file_counts

    # Iterate through chromosomes for this resolution
    for chrom in "${chromosomes[@]}"; do
        # Check if directory exists, if not, create it
        mkdir -p "$outdir/chr$chrom"

        # Initialize file creation counter for this chromosome at this resolution
        file_creation_count[$chrom]=0

        while IFS=$'\t' read -r idx base_name color; do
            expected_file="$outdir/chr$chrom/${base_name}_chr$chrom.txt"

            # Check if the expected file exists
            if [[ ! -f "$expected_file" ]]; then
                # File does not exist, so create it with default content
                echo -e "1\t1\t1" > "$expected_file"
                ((file_creation_count[$chrom]++))
                echo "$expected_file was created."
            fi
        done < "$filtered_list"
    done

    # Print the number of files created for each chromosome at this resolution
    echo "Files created for resolution $label:"
    for chrom in "${!file_creation_count[@]}"; do
        echo "Chromosome $chrom: ${file_creation_count[$chrom]}"
    done

    # Count the files in each chromosome directory and store in file_counts
    for chrom in "${chromosomes[@]}"; do
        file_counts[$chrom]=$(find "$outdir/chr$chrom" -type f | wc -l)
    done

    # Report the file counts for each chromosome directory at this resolution
    echo "Total file counts per chromosome directory for resolution $label:"
    for chrom in "${!file_counts[@]}"; do
        echo "Chromosome $chrom: ${file_counts[$chrom]}"
    done
done

# The check for equality in file counts across chromosomes is less meaningful when considering multiple resolutions,
# as the operation and file count are resolution-dependent.

# Convert resolutions string into an associative array for easy handling
declare -A resolutionMap
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    IFS=':' read -r -a keyValue <<< "$pair"
    resolutionMap[${keyValue[0]}]=${keyValue[1]}
done

# Base directory (assumed to be the main directory before creating subdirectories)
baseDir=$(pwd)

# Iterate through the resolutions and create directories
for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    imputeDir="hicluster_${label}_impute_dir"

    # Create the imputation directory for the current resolution
    mkdir -p "$imputeDir"
    echo "Created directory: $imputeDir"

    # Change into the imputation directory
    cd "$imputeDir"

    # Create directories for each chromosome within the imputation directory
    for i in {1..22}; do
        mkdir "chr$i"
        echo "Created directory: chr$i within $imputeDir"
    done

    # Change back to the main directory before processing the next resolution
    cd "$baseDir"
done

ls sc*.hic_matrix.txt.gz | perl -sne '
  use strict;
  use warnings;

  # Initialize variables with command line arguments
  our $software_directory;
  our $submission_script_directory;
  our $chrom_file;
  our $schicluster_env;
  our $resolutions_string;

  # Convert resolutions string to a hash where key is the resolution and value is the label
  my %resolutions = map { split /:/ } split /,/, $resolutions_string;

  # Remove newline and store the input filename
  chomp(my $in = $_);

  # Extract base ID by removing the file extension
  my $id = $in;
  $id =~ s/\.hic_matrix\.txt\.gz$//;

  # List of chromosomes to process
  my @chromosomes = (1..22);

  # Iterate through each resolution
  foreach my $res (keys %resolutions) {
      my $label = $resolutions{$res};

      # Iterate through each chromosome
      foreach my $chrom (@chromosomes) {
          # Construct the command for Hi-C imputation
          my $cmd = "source activate $schicluster_env && hicluster impute-cell ".
                    "--indir hicluster_${label}_raw_dir/chr${chrom}/ ".
                    "--outdir hicluster_${label}_impute_dir/chr${chrom}/ ".
                    "--cell $id ".
                    "--chrom chr$chrom ".
                    "--res $res ".
                    "--chrom_file $chrom_file";

          # Construct a job name
          my $job_name = "${label}_impute_${id}.${chrom}";
          $job_name =~ s/[^\w]/_/g;

          # Submit the job using the submission script directory path
          my $submit_cmd = "perl $submission_script_directory/ssubexc.pl ".
                           "\"$cmd\" \"$job_name\" 10 5";
          system($submit_cmd);
      }
  }
' -- -software_directory="$software_directory" -submission_script_directory="$submission_script_directory" -chrom_file="$chrom_file" -schicluster_env="$schicluster_env" -resolutions_string="$resolutions"

#########################################################################################                               ###compartment calling                                                                                                  #########################################################################################

cat > process_chromsizes.py <<EOF
import numpy as np
import pandas as pd

# Passed resolutions string from the shell script
resolutions_string = "$resolutions"

# Parse the resolutions string to a dictionary
resolutions = dict(res.split(":") for res in resolutions_string.split(","))

# Load chromosome sizes from a file specified by the chrom_file variable
chromsize = pd.read_csv('$chrom_file', sep='\t', header=None, index_col=0).to_dict()[1]

# Loop over each resolution and chromosome to create BED segments
for res_str, label in resolutions.items():
    res = int(res_str)  # Convert resolution string to integer
    for c in chromsize:
        # Calculate the number of segments based on chromosome size and resolution
        ngene = int(chromsize[c] // res) + 1
        # Generate BED format data for each segment
        bed = [[c, i*res, (i+1)*res] for i in range(ngene)]
        # Ensure the last segment extends to the end of the chromosome
        bed[-1][-1] = chromsize[c]
        # Save the BED data to a file, one for each chromosome
        outdir = f'./hicluster_{label}_impute_dir/bins'
        # Make sure the output directory exists
        import os
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        np.savetxt(f'{outdir}/{c}.bed', bed, fmt='%s', delimiter='\t')
        print(f'Resolution: {res_str}, Chromosome: {c}')
EOF

# Split the resolutions string into an array and iterate over it
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
    # Split each pair into resolution and label
    IFS=':' read -r resolution label <<< "$pair"
    # Create the output directory for the current resolution
    mkdir -p "./hicluster_${label}_impute_dir/bins"
    echo "Created directory: ./hicluster_${label}_impute_dir/bins"
done

# Execute the Python script                                                                                            
python process_chromsizes.py

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
    bedtools nuc -fi $software_directory/genomes/ucsc_hg19/hg19.fa -bed ./hicluster_${label}_impute_dir/bins/hg19.${label}_bin.bed -pattern CG -C > ./hicluster_${label}_impute_dir/bins/hg19.${label}_bin.cg_density.bed
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


#########################################################################################                               ###concatenate                                                                                                          #########################################################################################



