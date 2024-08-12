
######################################################################################                   #This shell script first checks to make sure all the previous preprocessing is complete.
#Then is runs the generate-matrix command from scHiCluster for every cell.
#If imputation boolean is set to True, then the imputed matrices from scHICluster 
#are computed as well. If the matrices are already present in the appropriate 
#directoryies, the computations are skipped. 
#Author: David Katz (davidkatz02@gmail.com)                                                              #####################################################################################

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
  hic_txt="$output_directory/${identifier}.hic_matrix.txt.gz"

  if [[ -f "$hic_txt" ]] && [[ -s "$hic_txt" ]]; then
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

# Parse the resolutions list into a hash
declare -A resolutionMap
for res in "${resolutions[@]}"; do
  IFS=':' read -r -a keyValue <<< "$res"
  resolutionMap[${keyValue[0]}]=${keyValue[1]}
done

# Process each Hi-C text file to generate matrix-cell
for hic_file in $output_directory/sc*.hic_matrix.txt.gz; do
  prefix=$(basename "$hic_file" .hic_matrix.txt.gz)

  for res in "${!resolutionMap[@]}"; do
    label=${resolutionMap[$res]}
    outdir="$output_directory/hic_${label}_raw_dir/"
    mkdir -p "$outdir"

    skip_computation=true
    for chrom in {1..22}; do
      output_matrix="$outdir/chr$chrom/${prefix}_chr${chrom}.txt"
      if [[ ! -f "$output_matrix" ]] || [[ ! -s "$output_matrix" ]]; then
        skip_computation=false
        break
      fi
    done

    if [ "$skip_computation" = true ]; then
      echo "Matrix for $prefix at resolution $label already exists for all chromosomes. Skipping computation."
      continue
    fi

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

if [ "$impute" = "True" ]; then
  echo "Imputation is enabled. Starting imputation steps..."
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
    echo "Processing file: $hic_file (prefix: $prefix)"

    for res in "${!resolutionMap[@]}"; do
      label=${resolutionMap[$res]}
      imputeDir="$output_directory/hicluster_${label}_impute_dir"
      echo "Resolution: $res, Label: $label, ImputeDir: $imputeDir"

      for chrom in {1..22}; do
        output_impute="$imputeDir/chr$chrom/${prefix}_chr${chrom}_pad1_std1_rp0.5_sqrtvc.hdf5"
        echo "Checking existence of imputed matrix: $output_impute"
        
        if [[ -f "$output_impute" ]] && [[ -s "$output_impute" ]]; then
          echo "Imputed matrix for $prefix for chromosome $chrom at resolution $label already exists. Skipping computation."
          continue
        fi

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
fi







