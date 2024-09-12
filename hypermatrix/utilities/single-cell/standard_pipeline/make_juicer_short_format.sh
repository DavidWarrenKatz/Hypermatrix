#make_juicer_short_format.sh

######################################################################################                   #This shell script first checks to make sure all the previous preprocessing is complete.
#Then is runs the generate-matrix command from scHiCluster for every cell.
#If imputation boolean is set to True, then the imputed matrices from scHICluster 
#are computed as well. If the matrices are already present in the appropriate 
#directoryies, the computations are skipped. 
#Author: David Katz (davidkatz02@gmail.com)                                                              #####################################################################################

#!/bin/bash

# Get the directory where this script (single_cell_pipeline.sh) is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Import the parameters from config.py (relative to the script's directory)
eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

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




