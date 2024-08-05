# Source the conda environment setup script
source /software/miniconda3/4.12.0/etc/profile.d/conda.sh

# Define environment names
schicluster_env=schicluster2
bisulfite_env=bisulfitehic27

# Activate the desired conda environment
conda activate $bisulfite_env

# Define directories and files
bam_directory='/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
software_directory='../../bin/softwarefiles'
chrom_file="../../bin/softwarefiles/hg19.autosome.chrom.sizes"
filtered_list='filtered_bam_list.txt'
output_directory='../../projects/single_cell_files'

# Create output directory if it doesn't exist
mkdir -p $output_directory

# Quality criteria
min_high_quality_reads=250000  

# Load samtools module
module load samtools

# Create or clear the filtered list file
> $filtered_list

# Count the initial number of BAM files
initial_bam_count=$(ls $bam_directory/sc*.b37.calmd.bam | wc -l)
echo "Initial number of BAM files: $initial_bam_count"

# Loop through BAM files and filter based on quality
for bam_file in $bam_directory/sc*.b37.calmd.bam; do
  echo "Evaluating $bam_file"
  
  # Count high-quality reads
  high_quality_reads=$(samtools view -c -q 30 -f 1 -F 1804 "$bam_file")
  
  # Check if the BAM file meets the quality criteria
  if (( high_quality_reads >= min_high_quality_reads )); then
    # Extract the identifier and add it to the filtered list
    identifier=$(basename "$bam_file" .b37.calmd.bam)
    echo "$identifier" >> $filtered_list
  fi
done

# Count the number of BAM files in the filtered list
filtered_bam_count=$(wc -l < $filtered_list)
echo "Number of BAM files in the filtered list: $filtered_bam_count"
echo "Filtered list of BAM files has been created."

# Create symbolic links to BAM files in the output directory based on the filtered list
while read identifier; do
  ln -s "$bam_directory/$identifier.b37.calmd.bam" "$output_directory/$identifier.bam"
done < $filtered_list

# Ensure that bisulfite_env_python_path is set correctly
bisulfite_env_python_path=$(which python)

# Process each BAM file in the output directory
for bam_file in $output_directory/sc*.bam; do
  echo "Processing $bam_file"
  
  # Extract prefix from filename for job naming and command construction
  prefix=$(basename "$bam_file" .bam)

  # Construct the command to be executed
  good_reads_bam="$output_directory/${prefix}.good_reads.bam"
  hic_txt="$output_directory/${prefix}.hic.txt"

  # Filter BAM file and convert to juicer_short format
  samtools view --threads 5 -bh -q 30 -f 1 -F 1804 "$bam_file" > "$good_reads_bam"
  if [ $? -ne 0 ]; then
    echo "Error processing $bam_file with samtools."
    continue
  fi

  "$bisulfite_env_python_path" "$software_directory/sam2juicer.py" -s "$good_reads_bam" -f "$bam_directory/hg19_DpnII.txt" > "$hic_txt"
  if [ $? -ne 0 ]; then
    echo "Error processing $good_reads_bam with sam2juicer.py."
    continue
  fi
  
  echo "Processed $bam_file -> $hic_txt"
done

echo "All BAM files have been processed."

