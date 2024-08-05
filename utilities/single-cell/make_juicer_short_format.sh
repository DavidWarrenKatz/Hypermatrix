#!/bin/bash

eval "$(python3 config_and_print.py)"

# Source the conda environment setup script
source /software/miniconda3/4.12.0/etc/profile.d/conda.sh

# Define environment names
schicluster_env=schicluster2
bisulfite_env=bisulfitehic27

# Activate the desired conda environment
conda activate $schicluster_env

#########################################################################################
### Generate contact matrix file for each single-cell from BAM file
### This block is resolution-independent; resolutions are specified in the next block
########################################################################################

# Check for successful file processing
successful_count=0
unsuccessful_files=()

while read identifier; do
  good_reads_bam="$output_directory/${identifier}.good_reads.bam"
  hic_txt="$output_directory/${identifier}.hic.txt"

  if [[ -f "$good_reads_bam" ]] && [[ -s "$good_reads_bam" ]] && [[ -f "$hic_txt" ]] && [[ -s "$hic_txt" ]]; then
    ((successful_count++))
  else
    unsuccessful_files+=("$identifier")
  fi
done < "$filtered_list"

echo "Number of successfully processed BAM files: $successful_count"
if [ ${#unsuccessful_files[@]} -eq 0 ]; then
  echo "All files were processed successfully."
else
  echo "Files not processed successfully:"
  for file in "${unsuccessful_files[@]}"; do
    echo "$file"
  done
fi

#########################################################################################
### Generate matrix-cell
### run the command from hic_cluster, since I might use their imputation step 
#########################################################################################

# This is a string of resolutions and names for the directory of each resolution
resolutions="1000000:1Mb"

echo "$resolutions"

# Parse the resolutions string into a hash
declare -A resolutionMap
IFS=',' read -r -a resolutionArray <<< "$resolutions"
for pair in "${resolutionArray[@]}"; do
  IFS=':' read -r -a keyValue <<< "$pair"
  resolutionMap[${keyValue[0]}]=${keyValue[1]}
done

# Process each Hi-C text file to generate matrix-cell
for hic_file in $output_directory/sc*.hic_matrix.txt.gz; do
  prefix=$(basename "$hic_file" .hic_matrix.txt.gz)

  for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    outdir="$output_directory/hicluster_${label}_raw_dir/"
    mkdir -p "$outdir"

    # Generate contact matrices for each resolution
    echo "Generating contact matrices for $prefix at resolution $label"

    hicluster generatematrix-cell \
      --infile "$hic_file" \
      --outdir "$outdir" \
      --chrom_file "$chrom_file" \
      --res "$res" \
      --cell "$prefix" \
      --chr1 1 \
      --pos1 2 \
      --chr2 5 \
      --pos2 6

    if [ $? -ne 0 ]; then
      echo "ERROR generating matrix for $prefix at resolution $label"
      continue
    fi

    echo "Generated matrix for $prefix at resolution $label"
  done
done

echo "All matrices have been generated."

#########################################################################################
### Impute-cell step
#########################################################################################

# Create directories for imputation results
for res in "${!resolutionMap[@]}"; do
  label=${resolutionMap[$res]}
  imputeDir="$output_directory/hicluster_${label}_impute_dir"
  mkdir -p "$imputeDir"
  for chrom in {1..22}; do
    mkdir -p "$imputeDir/chr$chrom"
  done
done

# Process each Hi-C matrix to perform imputation
for hic_file in $output_directory/sc*.hic_matrix.txt.gz; do
  prefix=$(basename "$hic_file" .hic_matrix.txt.gz)

  for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    imputeDir="$output_directory/hicluster_${label}_impute_dir"

    for chrom in {1..22}; do
      echo "Imputing $prefix for chromosome $chrom at resolution $label"

      hicluster impute-cell \
        --indir "$output_directory/hicluster_${label}_raw_dir/chr$chrom/" \
        --outdir "$imputeDir/chr$chrom/" \
        --cell "$prefix" \
        --chrom "chr$chrom" \
        --res "$res" \
        --chrom_file "$chrom_file"

      if [ $? -ne 0 ]; then
        echo "ERROR imputing $prefix for chromosome $chrom at resolution $label"
        continue
      fi

      echo "Imputed $prefix for chromosome $chrom at resolution $label"
    done
  done
done

echo "All imputation steps have been completed."

