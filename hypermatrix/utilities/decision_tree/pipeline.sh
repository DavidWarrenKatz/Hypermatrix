#!/bin/bash

# Define directories
input_dir="input_data_directory"
output_dir="output_data_directory"
trimmed_dir="${output_dir}/trimmed"
aligned_dir="${output_dir}/aligned"
counts_dir="${output_dir}/counts"

# Create output directories if they don't exist
mkdir -p "$trimmed_dir"
mkdir -p "$aligned_dir"
mkdir -p "$counts_dir"

# bcl2fastq command (you need to specify your own input and output parameters)
# bcl2fastq --input-dir $input_dir --output-dir $output_dir ...

# Trimmomatic for read trimming
for fq in ${input_dir}/*.fastq
do
    base=$(basename $fq .fastq)
    trimmomatic SE -phred33 \
        $fq \
        ${trimmed_dir}/${base}_trimmed.fastq \
        TRAILING:30 \
        MINLEN:20
done

# TopHat for alignment
# Replace with the path to your reference genome and annotation files
genome_index="path_to_your_genome_index"
annotation_file="path_to_your_annotation_file"

for trimmed_fq in ${trimmed_dir}/*.fastq
do
    base=$(basename $trimmed_fq _trimmed.fastq)
    tophat -o ${aligned_dir}/${base} \
        --no-novel-juncs \
        --read-mismatches 2 \
        --read-edit-dist 2 \
        --max-multihits 5 \
        $genome_index \
        $trimmed_fq
done

# Htseq-count for read counting
for bam in ${aligned_dir}/*.bam
do
    base=$(basename $bam .bam)
    htseq-count \
        -f bam \
        -s reverse \
        -t exon \
        -i gene_id \
        $bam \
        $annotation_file > ${counts_dir}/${base}_counts.txt
done

echo "RNA-Seq pipeline completed."

