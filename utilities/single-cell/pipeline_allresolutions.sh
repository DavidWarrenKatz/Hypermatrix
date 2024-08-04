######################################################################################
#This shell script performs the scHiCluster pipeline in single cell files
#For various resolutions
#Author: David Katz (davidkatz02@gmail.com) 
#####################################################################################

#This script assumes the hic files are in the form sc10.ACTTGA.b37.calmd.bam
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

bam_directory='/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
filtered_list='/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam/methylation/filter_low_qual/tsne/name.order.with_color.txt'
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
my $error_log="/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/new_processing9/logs/${random}ssubexc.err";
my $output_log="/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/new_processing9/logs/${random}ssubexc.out";

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
###this block is resolution independent, resolutions are specified in the next block
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

<<comment
#########################################################################################
###hicluster generatematrix-cell
###adjusts the resolution parameters as needed here, following the default format
#########################################################################################

#This is a string of resolutions and names fot the directory of each resolution
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

comment

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


#########################################################################################                               ###concatenate                                                                                                          #########################################################################################

#Here I am going to write the single chromosome concatenation

# concatenate cells for each chromosome
perl -e '@cs=(1..22);foreach $c(@cs){print "$c\n";`ls hicluster_250kb_impute_dir/chr${c}/*chr${c}_pad1_std1_rp0.5_sqrtvc.hdf5 > hicluster_100kb_chr${c}_impute_file_list.txt && hicluster embed-concatcell-chr --cell_list hicluster_250kb_chr${c}_impute_file_list.txt --outprefix hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr${c} --res 250000`;}'

perl -e '@cs=(2,12);foreach $c(@cs){`ls hicluster_250kb_impute_dir/chr${c}/*chr${c}_pad1_std1_rp0.5_sqrtvc.hdf5 > hicluster_250kb_chr${c}_impute_file_list.txt && hicluster embed-concatcell-chr --cell_list hicluster_250kb_chr${c}_impute_file_list.txt --outprefix hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr${c} --res 250000`;}'

ls ${impute_dir}/${chromosome}/*${imputation_mode}_${chromosome}.hdf5 > ${impute_file_list}
hicluster embed-concatcell-chr --cell_list ${impute_file_list} --outprefix ${embed_dir}/${imputation_mode}_${chromosome} --res ${resolution}



perl -e '@cs=(1..22);foreach $c(@cs){print "$c\n";`ls hicluster_250kb_impute_dir/chr${c}/*chr${c}_pad1_std1_rp0.5_sqrtvc.hdf5 > hicluster_250kb_chr${c}_impute_file_list.txt && hicluster embed-concatcell-chr --cell_list hicluster_250kb_chr${c}_impute_file_list.txt --outprefix hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_chr${c} --res 250000`;}'

# merge chromosomes (hicluster_100kb_embed_dir/pad1_std1_rp0.5_sqrtvc.svd50.hdf5)
ls hicluster_250kb_embed_dir/pad1_std1_rp0.5_sqrtvc_*npy > hicluster_250kb_embed_file_list.txt
hicluster embed-mergechr --embed_list hicluster_250kb_embed_file_list.txt --outprefix hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc

cat hicluster_250kb_chr1_impute_file_list.txt | perl -ne 'chomp;$in=$_;$in=~s/\S+\/(scD\w+\.\w+)_chr\S+.hdf5/$1/;print "$in\n";' > hicluster_250kb_sample_order.txt


######
# R code for last embedding step
#####

setwd("/jet/home/dnaase/startup/projects/hexa_seq/scNOMeHiC_RNA_IMR90_GM_20210701/bam/DNA/3D_genome")
library(Rtsne)
library(umap)
mc<-read.table("/jet/home/dnaase/startup/projects/hexa_seq/scNOMeHiC_RNA_IMR90_GM_20210701/bam/DNA/methy/filter_low_qual/tsne/name.order.HCG_methy.with_color.txt",sep="\t",header=F)
m_colors<-mc[,3]
names(m_colors)<-mc[,2]



#250kb (best, with top 3 SVD)
ind<-read.table("hicluster_250kb_sample_order.txt",sep="\t",header=F)
mat<-read.table("/jet/home/dnaase/startup/projects/hexa_seq/scNOMeHiC_RNA_IMR90_GM_20210701/bam/DNA/3D_genome/hicluster_250kb_embed_dir/all_merged.pad1_std1_rp0.5_sqrtvc.svd50.txt",sep="\t",header=F)
rownames(mat)<-ind[,1]

set.seed(42) # Sets seed for reproducibility
for(j in c(5,10,15,20,25,30)){
	pdf(paste("hicluster_3d_svd50.250kb.plexity_",j,".pdf",sep=""),paper="special", height=5, width=5)
	tsne_out <- Rtsne(as.matrix(mat[,1:3]),perplexity = j, max_iter=5000) # Run TSNE
	plot(tsne_out$Y,col=m_colors[rownames(mat)],pch=16) # Plot the result
	dev.off()
}

pdf("hicluster_3d_svd50.250kb.umap.pdf",paper="special", height=5, width=5)
umap_out <- umap(as.matrix(mat[,1:3])) # Run umap
plot(umap_out$layout,col=m_colors[rownames(mat)],pch=16) # Plot the result
dev.off()





