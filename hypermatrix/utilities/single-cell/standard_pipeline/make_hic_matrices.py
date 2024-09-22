#!/usr/bin/env python3
# File: make_hic_matrices.py

import argparse
import gzip
import h5py
import json
import logging
import numpy as np
import os
import sys
import glob
from scipy.sparse import csr_matrix, triu, tril
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

def load_hic_data(filepath):
    """Load Hi-C data from a text file."""
    data = []
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                i, j, reads = int(parts[0]), int(parts[1]), float(parts[2])
                data.append((i, j, reads))
    return data

def create_matrix(data, size):
    """Creates a symmetric matrix from Hi-C data."""
    if not data:
        return None
    matrix = csr_matrix((size, size), dtype=float)
    rows = []
    cols = []
    values = []
    for i, j, reads in data:
        rows.append(i)
        cols.append(j)
        values.append(reads)
        if i != j:
            rows.append(j)
            cols.append(i)
            values.append(reads)
    matrix = csr_matrix((values, (rows, cols)), shape=(size, size))
    return matrix

def emphasize_interactions(matrix, max_distance):
    """Emphasize interactions by adding offsets of contacts."""
    emphasized_matrix = matrix.copy()
    for offset in range(1, max_distance + 1):
        emphasized_matrix += triu(matrix, k=offset) + tril(matrix, k=-offset)
    return emphasized_matrix

def process_matrices(input_dir, output_emphasized_dir, max_distance, size):
    """Process each Hi-C data file, compute matrices, and save the results."""
    os.makedirs(output_emphasized_dir, exist_ok=True)
    logging.info(f"Processing input directory: {input_dir}")
    for file_path in glob.glob(os.path.join(input_dir, '*.txt')):
        logging.info(f"Processing file: {file_path}")

        file_name = os.path.splitext(os.path.basename(file_path))[0] + '.h5'
        output_emphasized_path = os.path.join(output_emphasized_dir, file_name)

        # Skip computation if the output file already exists
        if os.path.exists(output_emphasized_path):
            logging.info(f"Skipping {output_emphasized_path}, emphasized Hi-C matrix already exists.")
            continue

        data = load_hic_data(file_path)
        if not data:
            logging.warning(f"No data loaded from {file_path}")
            continue
        csr_mat = create_matrix(data, size)
        if csr_mat is None:
            logging.warning(f"Failed to create matrix from data in {file_path}")
            continue

        emphasized_matrix = emphasize_interactions(csr_mat, max_distance)

        with h5py.File(output_emphasized_path, 'w') as output_file:
            grp = output_file.create_group('Matrix')
            grp.create_dataset('data', data=emphasized_matrix.data)
            grp.create_dataset('indices', data=emphasized_matrix.indices)
            grp.create_dataset('indptr', data=emphasized_matrix.indptr)
            grp.attrs['shape'] = emphasized_matrix.shape

        logging.info(f"Saved emphasized Hi-C matrix to {output_emphasized_path}")

def main():
    parser = argparse.ArgumentParser(description='Generate Hi-C matrices and tensors.')
    parser.add_argument('--config', type=str, default='config.json', help='Path to the configuration JSON file')
    parser.add_argument('--bam_directory', type=str, help='Path to the Hi-C BAM directory')
    parser.add_argument('--chrom_file', type=str, help='Path to the chromosome sizes file')
    parser.add_argument('--resolutions', type=str, help='Resolution in the format value:label, e.g., "1000000:1Mb"')
    parser.add_argument('--output_directory', type=str, help='Output directory')
    parser.add_argument('--chromosomes', type=str, nargs='+', help='List of chromosomes to process')
    args = parser.parse_args()

    # Load config
    config_path = args.config
    config = load_config(config_path)

    # Override config with command-line arguments if provided
    for key in ['bam_directory', 'chrom_file', 'resolutions', 'output_directory', 'chromosomes']:
        arg_value = getattr(args, key)
        if arg_value is not None:
            config[key] = arg_value

    # Access values from config
    bam_directory = config.get('bam_directory')
    chrom_file = config.get('chrom_file')
    resolutions = config.get('resolutions')
    output_directory = config.get('output_directory')
    chromosomes = config.get('chromosomes')

    # Validate required parameters
    required_params = {
        'bam_directory': bam_directory,
        'chrom_file': chrom_file,
        'resolutions': resolutions,
        'output_directory': output_directory,
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

    # Process resolutions
    for resolution_str in resolutions:
        try:
            resolution, resolution_label = parse_resolution(resolution_str)
        except ValueError as e:
            logging.error(e)
            sys.exit(1)

        max_genomic_distance = int(10_000_000 / resolution) + 1
        base_input_dir = os.path.join(output_directory, f'hic_{resolution_label}_raw_dir')
        base_output_emphasized_dir = os.path.join(output_directory, f'hic_{resolution_label}_emphasized_dir')

        for chrom in chromosomes:
            input_dir = os.path.join(base_input_dir, chrom)
            output_emphasized_dir = os.path.join(base_output_emphasized_dir, chrom)

            # Calculate the number of bins for the chromosome
            chrom_length = chromosome_lengths.get(chrom)
            if chrom_length is None:
                logging.warning(f"Chromosome {chrom} not found in chromosome sizes file.")
                continue
            size = chrom_length // resolution + 1  # Number of bins

            logging.info(f'Processing chromosome {chrom} at resolution {resolution_label}')
            process_matrices(input_dir, output_emphasized_dir, max_genomic_distance, size)

    logging.info("Processing completed successfully.")

if __name__ == '__main__':
    main()
