#!/usr/bin/env python3
# File: make_methy_matrices.py

import argparse
import gzip
import h5py
import json
import logging
import numpy as np
import os
import sys
from collections import defaultdict

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config(config_path):
    if not os.path.exists(config_path):
        logging.error(f"Configuration file not found at {config_path}")
        sys.exit(1)
    with open(config_path, 'r') as f:
        return json.load(f)

def parse_resolution(resolution_str):
    if ':' in resolution_str:
        resolution_value, resolution_label = resolution_str.split(':')
        try:
            resolution = int(resolution_value)
            return resolution, resolution_label
        except ValueError:
            raise ValueError(f"Resolution value should be an integer: '{resolution_value}' in '{resolution_str}'")
    else:
        raise ValueError(f"Invalid resolution format: '{resolution_str}'. Expected format 'value:label', e.g., '1000000:1Mb'.")

def calculate_matrix(file_path, chromosomes, resolution):
    # Initialize a dictionary to hold binned methylation values per chromosome
    binned_data = defaultdict(lambda: defaultdict(list))  # {chrom: {bin_index: [methylation_levels]}}

    with gzip.open(file_path, 'rt') as file:
        for line in file:
            fields = line.strip().split('\t')
            # Assuming columns: chrom, start, end, methylation_level
            if len(fields) < 4:
                continue
            chrom = fields[0]
            if chrom not in chromosomes:
                continue  # Skip chromosomes we're not interested in
            try:
                start = int(fields[1])
                methylation_level = float(fields[3])
            except ValueError:
                continue  # Skip lines with invalid data
            bin_index = start // resolution
            binned_data[chrom][bin_index].append(methylation_level)

    # Calculate mean methylation levels per bin
    methylation_matrices = {}
    for chrom in binned_data:
        bins = sorted(binned_data[chrom].keys())
        methylation_levels = []
        for bin_index in bins:
            levels = binned_data[chrom][bin_index]
            mean_level = np.mean(levels)
            methylation_levels.append(mean_level)
        methylation_matrices[chrom] = np.array(methylation_levels)

    return methylation_matrices

def main():
    parser = argparse.ArgumentParser(description='Generate methylation matrices and tensors.')
    parser.add_argument('--config', type=str, default='config.json', help='Path to the configuration JSON file')
    parser.add_argument('--methy_directory', type=str, help='Path to the methylation directory')
    parser.add_argument('--chrom_file', type=str, help='Path to the chromosome sizes file')
    parser.add_argument('--resolutions', type=str, help='Resolution in the format value:label, e.g., "1000000:1Mb"')
    parser.add_argument('--output_directory', type=str, help='Output directory')
    parser.add_argument('--reference_genome', type=str, help='Reference genome (e.g., hg19 or hg38)')
    parser.add_argument('--methy_data_file', type=str, help='Name of the methylation data file')
    parser.add_argument('--chromosomes', type=str, nargs='+', help='List of chromosomes to process')
    args = parser.parse_args()

    # Load config
    config_path = args.config
    config = load_config(config_path)

    # Override config with command-line arguments if provided
    for key in ['methy_directory', 'chrom_file', 'resolutions', 'output_directory', 'reference_genome', 'methy_data_file', 'chromosomes']:
        arg_value = getattr(args, key)
        if arg_value is not None:
            config[key] = arg_value

    # Access values from config
    methy_directory = config.get('methy_directory')
    chrom_file = config.get('chrom_file')
    resolutions = config.get('resolutions')
    output_directory = config.get('output_directory')
    reference_genome = config.get('reference_genome')
    methy_data_file = config.get('methy_data_file')
    chromosomes = config.get('chromosomes')

    # Validate required parameters
    required_params = {
        'methy_directory': methy_directory,
        'chrom_file': chrom_file,
        'resolutions': resolutions,
        'output_directory': output_directory,
        'reference_genome': reference_genome,
        'methy_data_file': methy_data_file,
        'chromosomes': chromosomes
    }
    for param_name, param_value in required_params.items():
        if not param_value:
            logging.error(f"Missing required parameter: {param_name}")
            sys.exit(1)

    # Convert chromosomes to strings if they are integers
    chromosomes = [str(chrom) for chrom in chromosomes]

    # Log configuration
    logging.info(f"Using configuration: {json.dumps(required_params, indent=2)}")

    # Ensure resolutions is a list
    if isinstance(resolutions, str):
        resolutions = [resolutions]

    # Process resolutions
    for resolution_str in resolutions:
        try:
            resolution, resolution_label = parse_resolution(resolution_str)
        except ValueError as e:
            logging.error(e)
            sys.exit(1)

        # Update paths
        methy_output_dir = os.path.join(output_directory, f"methy_{resolution_label}_outerproduct_dir")
        methy_matrix_path = os.path.join(methy_directory, methy_data_file)

        # Ensure output directory exists
        os.makedirs(methy_output_dir, exist_ok=True)

        # Read chromosome sizes
        chromosome_lengths = {}
        try:
            with open(chrom_file, 'r') as f:
                for line in f:
                    chrom, length = line.strip().split()
                    if chrom in chromosomes:
                        chromosome_lengths[chrom] = int(length)
        except Exception as e:
            logging.error(f"Error reading chromosome sizes: {e}")
            sys.exit(1)

        # Calculate methylation matrices per chromosome
        try:
            logging.info(f"Calculating methylation matrices from {methy_matrix_path}")
            methylation_matrices = calculate_matrix(methy_matrix_path, chromosomes, resolution)
        except Exception as e:
            logging.error(f"Error calculating methylation matrices: {e}")
            sys.exit(1)

        # Process each chromosome
        for chrom in chromosomes:
            if chrom not in methylation_matrices:
                logging.warning(f"No data found for chromosome {chrom}")
                continue

            data_vector = methylation_matrices[chrom]
            logging.info(f'Chromosome {chrom}: Data vector length: {len(data_vector)}')

            # Center the data
            vector_centered = data_vector - np.mean(data_vector)

            # Compute outer product (if size is manageable)
            vector_length = len(vector_centered)
            if vector_length > 10000:
                logging.warning(f"Chromosome {chrom} has {vector_length} bins, which may cause memory issues.")
                logging.info(f"Skipping outer product computation for chromosome {chrom} due to size.")
                continue

            outer_product_matrix = np.outer(vector_centered, vector_centered)

            # Save individual matrix
            chromosome_dir = os.path.join(methy_output_dir, chrom)
            os.makedirs(chromosome_dir, exist_ok=True)
            individual_matrix_filename = os.path.join(chromosome_dir, f"{chrom}_outer_product.h5")
            with h5py.File(individual_matrix_filename, 'w') as individual_matrix_file:
                individual_matrix_file.create_dataset('methylation', data=outer_product_matrix, compression="gzip")

            logging.info(f"Saved outer product matrix for chromosome {chrom}")

        logging.info(f"Processing for resolution {resolution_label} completed.")

    logging.info("Processing completed successfully.")

if __name__ == '__main__':
    main()
