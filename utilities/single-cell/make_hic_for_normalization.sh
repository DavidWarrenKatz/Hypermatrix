#!/bin/bash

# Configuration variables
eval "$(python3 config_and_print.py)"
juicer_tools_path="$software_directory/juicer_tools_1.22.01.jar"

# Extract resolution and label from the resolutions string
IFS=':' read -r resolution resolution_label <<< "$resolutions"

data_path="$output_directory/hic_${resolution_label}_raw_dir"

# List of autosomal chromosomes
chromosomes=$(seq 1 22)  # Chromosomes 1 to 22

# Read prefixes from the filtered list file
datasets=()
while IFS= read -r line; do
    datasets+=("$line")
done < "$filtered_list"

# Function to convert to HiC Short format
convert_to_hic_short_format() {
    input_file=$1
    chrom=$2
    output_file=$3
    resolution=$4

    if [[ -f "$output_file" ]]; then
        echo "Formatted file already exists: $output_file. Skipping conversion."
        return
    fi

    while IFS= read -r line; do
        parts=($line)
        if [[ ${#parts[@]} -ne 3 ]]; then
            echo "Skipping invalid line in $input_file: $line"
            continue
        fi

        pos1=$((${parts[0]} * resolution))
        pos2=$((${parts[1]} * resolution))
        score=${parts[2]}

        echo "0 $chrom $pos1 0 0 $chrom $pos2 1 $score" >> "$output_file"
    done < "$input_file"

    echo "Conversion completed: $output_file"
}

# Function to create .hic file
create_hic_file() {
    chrom=$1
    dataset=$2

    short_form_file="$data_path/chr$chrom/${dataset}_chr${chrom}.txt"
    formatted_file="$data_path/chr$chrom/${dataset}_chr${chrom}_juicer_formatted.txt"
    hic_file="$data_path/chr$chrom/${dataset}_chr${chrom}.hic"

    # Convert to HiC Short format with score using the provided resolution
    convert_to_hic_short_format "$short_form_file" "$chrom" "$formatted_file" "$resolution"

    # Run juicer tools to create .hic file
    if ! validate_hic_file "$hic_file" "$chrom"; then
        java -Xmx100G -jar "$juicer_tools_path" pre -r "$resolution" "$formatted_file" "$hic_file" hg19
    fi

    # Retry mechanism
    max_retries=10
    for ((i=1; i<=max_retries; i++)); do
        if validate_hic_file "$hic_file" "$chrom"; then
            echo "HIC file successfully created: $hic_file"
            return 0
        else
            echo "HIC file not created yet or is invalid, retrying ($i/$max_retries)..."
            sleep 20
        fi
    done

    echo "Failed to create a valid HIC file after $max_retries retries: $hic_file"
    return 1
}

# Function to validate .hic file using juicer_tools
validate_hic_file() {
    hic_file=$1
    chrom=$2

    java -jar "$juicer_tools_path" dump observed NONE "$hic_file" "$chrom" "$chrom" BP "$resolution" /dev/null
    return $?
}

# Function to apply KR normalization and convert back to short format
kr_normalize_and_convert() {
    hic_file=$1
    chrom=$2
    dataset=$3

    if [[ ! -f "$hic_file" ]]; then
        echo "Skipping KR normalization due to missing or invalid HIC file for ${dataset} chr${chrom}."
        return
    fi

    kr_matrix_file="$data_path/chr$chrom/${dataset}_chr${chrom}_KR_matrix.txt"
    kr_short_file="$data_path/chr$chrom/${dataset}_chr${chrom}_KR_short.txt"

    java -jar "$juicer_tools_path" dump observed KR "$hic_file" "chr$chrom" "chr$chrom" BP "$resolution" "$kr_matrix_file"

    awk '{print $1 "\t" $2 "\t" $3}' "$kr_matrix_file" > "$kr_short_file"

    echo "KR normalization completed for $dataset chr${chrom}. Short format saved to $kr_short_file"
}

# Main execution flow

# Phase 1: Create all .hic files
for chrom in $chromosomes; do
    for dataset in "${datasets[@]}"; do
        create_hic_file "$chrom" "$dataset"
    done
done

# Phase 2: Apply KR normalization to all created .hic files
for chrom in $chromosomes; do
    for dataset in "${datasets[@]}"; do
        hic_file="$data_path/chr$chrom/${dataset}_chr${chrom}.hic"
        kr_normalize_and_convert "$hic_file" "$chrom" "$dataset"
    done
done

echo "Processing completed."

