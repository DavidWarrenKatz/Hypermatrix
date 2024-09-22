#!/usr/bin/env python3
# File: config.py

import os
import sys
import json
import logging
import scipy.io

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Root directory based on config.py location
root_dir = os.path.dirname(os.path.realpath(__file__))

# Path to the configuration JSON file
config_json_path = os.path.join(root_dir, 'config.json')

def load_config():
    """
    Load configuration from a JSON file.
    """
    if not os.path.exists(config_json_path):
        logging.error(f"Configuration file not found at {config_json_path}")
        sys.exit(1)
    with open(config_json_path, 'r') as f:
        return json.load(f)

# Load configuration
config = load_config()

# Validate required parameters
required_params = [
    'bam_directory',
    'methy_directory',
    'output_directory',
    'reference_genome',
    'genome_fa',    # Path to the genome .fa file
    'chrom_file'    # Path to the genome chrom.sizes file
]

missing_params = [param for param in required_params if not config.get(param)]
if missing_params:
    logging.error(f"Missing required configuration parameters: {', '.join(missing_params)}")
    sys.exit(1)

# Set up genome directories
genomes_directory = os.path.join(root_dir, 'genomes')
os.makedirs(genomes_directory, exist_ok=True)

# Reference genome and paths
reference_genome = config['reference_genome']
genome_fa = config['genome_fa']
chrom_file = config['chrom_file']

def setup_genome_files(genome, genome_fa, chrom_file):
    """
    Validate the existence of genome files.
    """
    if not os.path.exists(genome_fa):
        logging.error(f"Genome FASTA file not found: {genome_fa}")
        sys.exit(1)
    if not os.path.exists(chrom_file):
        logging.error(f"Chromosome sizes file not found: {chrom_file}")
        sys.exit(1)
    return genome_fa, chrom_file

# Validate genome files
genome_fa, chrom_file = setup_genome_files(reference_genome, genome_fa, chrom_file)

# Path for the restriction enzyme fragment file (assuming DpnII was used)
fragments_file = os.path.join(genomes_directory, f"{reference_genome}_DpnII.txt")

# Update config with derived parameters
config.update({
    "software_directory": genomes_directory,  # Pointing to the genomes directory
    "fragments_file": fragments_file
})

# Save the configuration dictionary for later use (e.g., in MATLAB)
config_mat_path = os.path.join(root_dir, 'config.mat')
scipy.io.savemat(config_mat_path, config)
