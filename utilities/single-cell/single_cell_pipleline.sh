source /software/miniconda3/4.12.0/etc/profile.d/conda.sh

schicluster_env=schicluster2
bisulfite_env=bisulfitehic27

conda activate $bisulfite_env

bam_directory='/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam'
software_directory='/home/dwk681/workspace/softwareFiles'
chrom_file="/home/dwk681/workspace/Genomic_Files/ucsc_hg19/hg19.autosome.chrom.sizes"
filtered_list='filtered_bam_list.txt'

# Quality criteria
min_high_quality_reads=100000  

# Load samtools module
module load samtools

# Create or clear the filtered list file
> $filtered_list

# Evaluate each BAM file for quality and generate the filtered list
for bam_file in $bam_directory/sc*.b37.calmd.bam; do
  # Count high-quality reads
  high_quality_reads=$(samtools view -c -q 30 -f 1 -F 1804 "$bam_file")
  
  # Check if the BAM file meets the quality criteria
  if (( high_quality_reads >= min_high_quality_reads )); then
    # Extract the identifier and add it to the filtered list
    identifier=$(basename "$bam_file" .b37.calmd.bam)
    echo "$identifier" >> $filtered_list
  fi
done

echo "Filtered list of BAM files has been created."

# Create symbolic links to BAM files in the current directory based on the filtered list
while read identifier; do
  ln -s "$bam_directory/$identifier.b37.calmd.bam" "./$identifier.bam"
done < $filtered_list

# Iterate over BAM files matching the pattern "sc*.bam"
for bam_file in sc*.bam; do
  # Extract prefix from filename for job naming and command construction
  prefix=$(basename "$bam_file" .bam)

  # Construct the command to be executed
  good_reads_bam="${prefix}.good_reads.bam"
  hic_txt="${prefix}.hic.txt"

  # Filter BAM file and convert to juicer_short format
  samtools view --threads 5 -bh -q 30 -f 1 -F 1804 "$bam_file" > "$good_reads_bam"
  "$bisulfite_env_python_path" "$software_directory/sam2juicer.py" -s "$good_reads_bam" -f "$bam_directory/hg19_DpnII.txt" > "$hic_txt"
  
  echo "Processed $bam_file -> $hic_txt"
done

echo "All BAM files have been processed."

