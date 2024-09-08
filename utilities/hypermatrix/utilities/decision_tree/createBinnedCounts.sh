#!/bin/bash

# Define directories
input_dir="input_data_directory"
output_dir="output_data_directory"
trimmed_dir="${output_dir}/trimmed"
aligned_dir="${output_dir}/aligned"
counts_dir="${output_dir}/counts"
binned_counts_dir="${output_dir}/binned_counts"

# Define bin size (n)
bin_size=YOUR_BIN_SIZE  # Replace YOUR_BIN_SIZE with the desired bin size

# Create output directories if they don't exist
mkdir -p "$trimmed_dir"
mkdir -p "$aligned_dir"
mkdir -p "$counts_dir"
mkdir -p "$binned_counts_dir"

# [Your existing pipeline steps for bcl2fastq, Trimmomatic, TopHat, and Htseq-count]

# Assuming TopHat outputs .bam files and a sorted index (.bai) exists for each
for bam in ${aligned_dir}/*.bam
do
    base=$(basename $bam .bam)

    # Generate binned counts
    bedtools genomecov -ibam $bam -dz -g path_to_your_genome_file > ${binned_counts_dir}/${base}_genomecov.txt

    # Divide the genome into bins and count reads per bin
    awk -v binSize=$bin_size '
    BEGIN {currentBin = 1; binCount = 0;}
    {
        binIndex = int(($2 + $3) / 2 / binSize) + 1;
        if (binIndex > currentBin) {
            print currentBin, binCount;
            currentBin = binIndex;
            binCount = 0;
        }
        binCount += $4;
    }
    END {print currentBin, binCount;}' ${binned_counts_dir}/${base}_genomecov.txt > ${binned_counts_dir}/${base}_binned_counts.txt
done

echo "RNA-Seq pipeline completed with binned counts."

