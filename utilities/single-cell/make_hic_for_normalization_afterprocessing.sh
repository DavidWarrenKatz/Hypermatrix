#!/bin/bash

# Configuration variables
eval "$(python3 config_and_print.py)"
juicer_tools_path="$software_directory/juicer_tools_1.22.01.jar"

# Log file
log_file="hic_processing.log"
exec > >(tee -i $log_file)
exec 2>&1

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

# Retry and sleep configuration
max_retries=3
sleep_time=5

# Function to convert to HiC Short format
convert_to_hic_short_format() {
    input_file=$1
    chrom=$2
    output_file=$3
    resolution=$4

    if [[ -f "$output_file" && -s "$output_file" ]]; then
        echo "Formatted file already exists and is non-empty: $output_file. Skipping conversion."
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

# Function to create .hic file with retry mechanism
create_hic_file() {
    chrom=$1
    dataset=$2

    short_form_file="$data_path/chr$chrom/${dataset}_chr${chrom}_processed.txt"
    formatted_file="$data_path/chr$chrom/${dataset}_chr${chrom}_processed_juicer_formatted.txt"
    hic_file="$data_path/chr$chrom/${dataset}_chr${chrom}_processed.hic"

    convert_to_hic_short_format "$short_form_file" "$chrom" "$formatted_file" "$resolution"

    for ((i=1; i<=max_retries; i++)); do
        if ! validate_hic_file "$hic_file" "$chrom"; then
            java -Xmx100G -jar "$juicer_tools_path" pre -r "$resolution" "$formatted_file" "$hic_file" hg19
        fi
        if validate_hic_file "$hic_file" "$chrom"; then
            echo "HIC file successfully created: $hic_file"
            return 0
        else
            echo "HIC file not created yet or is invalid, retrying ($i/$max_retries)..."
            sleep $((sleep_time * i))
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

# Function to perform KR normalization with retry mechanism
kr_normalize() {
    hic_file=$1
    chrom=$2
    dataset=$3

    kr_bp_file="$data_path/chr$chrom/${dataset}_chr${chrom}_KR_short_bp.txt"

    if [[ -f "$kr_bp_file" && -s "$kr_bp_file" ]]; then
        echo "KR BP matrix already exists and is non-empty: $kr_bp_file. Skipping KR normalization."
        return
    fi

    for ((i=1; i<=max_retries; i++)); do
        java -jar "$juicer_tools_path" dump observed KR "$hic_file" "chr$chrom" "chr$chrom" BP "$resolution" "$kr_bp_file"
        if [[ -f "$kr_bp_file" && -s "$kr_bp_file" ]]; then
            echo "KR normalization completed for $dataset chr${chrom}. Matrix saved to $kr_bp_file"
            return 0
        else
            echo "KR normalization not completed yet or is invalid, retrying ($i/$max_retries)..."
            sleep $((sleep_time * i))
        fi
    done

    echo "Failed to create a valid KR matrix after $max_retries retries: $kr_bp_file"
    return 1
}

# Function to convert KR normalized matrix to short format with retry mechanism
convert_kr_to_short_format() {
    kr_bp_file=$1
    chrom=$2
    dataset=$3
    resolution=$4

    kr_bin_file="$data_path/chr$chrom/${dataset}_chr${chrom}_KR_short_bin.txt"

    if [[ -f "$kr_bin_file" && -s "$kr_bin_file" ]]; then
        echo "KR short format file already exists and is non-empty: $kr_bin_file. Skipping conversion."
        return
    fi

    for ((i=1; i<=max_retries; i++)); do
        awk -v res="$resolution" '{print int($1/res) "\t" int($2/res) "\t" $3}' "$kr_bp_file" > "$kr_bin_file"
        if [[ -f "$kr_bin_file" && -s "$kr_bin_file" ]]; then
            echo "Conversion to short format completed for $dataset chr${chrom}. Short format saved to $kr_bin_file"
            return 0
        else
            echo "Conversion to short format not completed yet or is invalid, retrying ($i/$max_retries)..."
            sleep $((sleep_time * i))
        fi
    done

    echo "Failed to create a valid KR short format file after $max_retries retries: $kr_bin_file"
    return 1
}

# Main execution flow

# Phase 1: Create all .hic files
for chrom in $chromosomes; do
    for dataset in "${datasets[@]}"; do
        create_hic_file "$chrom" "$dataset" &
    done
    wait  # Wait for all background tasks to finish before proceeding
done

# Phase 2: Perform KR normalization on all .hic files
for chrom in $chromosomes; do
    for dataset in "${datasets[@]}"; do
        hic_file="$data_path/chr$chrom/${dataset}_chr${chrom}_processed.hic"
        kr_normalize "$hic_file" "$chrom" "$dataset" &
    done
    wait  # Wait for all background tasks to finish before proceeding
done

# Phase 3: Convert KR normalized matrices to short format
for chrom in $chromosomes; do
    for dataset in "${datasets[@]}"; do
        kr_bp_file="$data_path/chr$chrom/${dataset}_chr${chrom}_KR_short_bp.txt"
        convert_kr_to_short_format "$kr_bp_file" "$chrom" "$dataset" "$resolution" &
    done
    wait  # Wait for all background tasks to finish before proceeding
done

echo "Processing completed."

# Optional: Clean-up temporary files
# rm -f "$data_path"/chr*/${dataset}_chr*_{juicer_formatted,KR_matrix}.txt


















