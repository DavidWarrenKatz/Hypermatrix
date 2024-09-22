#!/usr/bin/env python3
# File: hypermatrix/main.py

import argparse
import datetime
import json
import logging
import os
import subprocess
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Set the script directory and config file path
script_dir = os.path.dirname(os.path.realpath(__file__))
config_filename = 'config.json'
config_path = os.path.join(script_dir, config_filename)
sys.path.append(script_dir)

def load_config():
    if not os.path.exists(config_path):
        logging.error(f"Configuration file not found at {config_path}")
        sys.exit(1)
    with open(config_path, 'r') as f:
        return json.load(f)

def save_config(config_data):
    with open(config_path, 'w') as f:
        json.dump(config_data, f, indent=4)
    logging.info(f"Configuration updated in {config_path}")

VERSION = """
Hypermatrix version 0.1 - A tool for integrating
multi-omics data and epigenetic analysis using
advanced tensor techniques.
"""

def main():
    global config
    config = load_config()

    parser = argparse.ArgumentParser(description='Hypermatrix command-line tool')
    parser.add_argument('-v', '--version', action='version', version=VERSION)

    # Subcommands
    subparsers = parser.add_subparsers(dest='command', required=True)

    # preprocess subcommand
    preprocess_parser = subparsers.add_parser('preprocess', help='Preprocess BAM files for analysis')
    group = preprocess_parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--nomehic', action='store_true', help='Process BAM files from the scNOMe-HiC technique')
    group.add_argument('--m3C', action='store_true', help='Process BAM files from the m3C-seq technique')
    preprocess_parser.add_argument('--input_dir', type=str, required=True, help='Directory containing BAM files')
    preprocess_parser.add_argument('--output_dir', type=str, required=True, help='Directory where the output files will be saved')
    preprocess_parser.add_argument('--ref_genome', type=str, required=True, help='Reference genome (e.g., hg19 or hg38)')
    preprocess_parser.add_argument('--genome_fa', type=str, required=True, help='Path to the genome .fa file')
    preprocess_parser.add_argument('--chrom_sizes', type=str, required=True, help='Path to the genome chrom.sizes file')

    # ABcluster subcommand
    abcluster_parser = subparsers.add_parser('ABcluster', help='Run the ABcluster analysis')
    abcluster_parser.add_argument('-m', '--methy', type=str, help='Path to the single-cell CpG methylation directory')
    abcluster_parser.add_argument('-i', '--hic', type=str, help='Path to the single-cell chromosome conformation directory')
    abcluster_parser.add_argument('-o', '--output_dir', type=str, help='Directory where the output files will be saved')
    abcluster_parser.add_argument('-c', '--cumulant', action='store_true', help='Run cumulant shell script')
    abcluster_parser.add_argument('-p', '--impute', action='store_true', help='Run impute shell script')
    abcluster_parser.add_argument('-g', '--genome_id', type=str, help='Reference genome (e.g., hg19 or hg38)')
    abcluster_parser.add_argument('-r', '--res', type=str, help='Resolution, size of genomic bins')
    abcluster_parser.add_argument('--genome_fa', type=str, required=True, help='Path to the genome .fa file')
    abcluster_parser.add_argument('--chrom_sizes', type=str, required=True, help='Path to the genome chrom.sizes file')

    # differentiate_chromosomes subcommand
    diffchrom_parser = subparsers.add_parser('differentiate_chromosomes', help='Analyze distinct A/B compartments for homologous chromosomes')
    diffchrom_parser.add_argument('--hic_file', type=str, help='Path to the Hi-C data file')
    diffchrom_parser.add_argument('--epigenetic_file', type=str, help='Path to the epigenetic data file (optional)')
    diffchrom_parser.add_argument('--output_dir', type=str, help='Directory where the output files will be saved')

    args = parser.parse_args()

    # Dispatch the command
    if args.command == 'ABcluster':
        abcluster(args)
    elif args.command == 'differentiate_chromosomes':
        diffchrom(args)
    elif args.command == 'preprocess':
        preprocess(args)
    else:
        logging.error(f"Unknown command: {args.command}. Try 'ABcluster', 'differentiate_chromosomes', or 'preprocess'.")
        parser.print_help()
        sys.exit(1)

def abcluster(args):
    global config
    logging.info("Running ABcluster module")

    # Collect updates for the config file
    updates = {}
    if args.methy:
        updates['methy_directory'] = args.methy
    if args.hic:
        updates['bam_directory'] = args.hic
    if args.output_dir:
        updates['output_directory'] = args.output_dir
    if args.genome_id:
        updates['reference_genome'] = args.genome_id
    if args.res:
        updates['resolutions'] = [args.res]
    if args.genome_fa:
        updates['genome_fa'] = args.genome_fa
    if args.chrom_sizes:
        updates['chrom_file'] = args.chrom_sizes

    # Apply the updates to the config
    config.update(updates)
    save_config(config)

    # Ensure output directory exists
    output_dir = config.get('output_directory')
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Decide which script to run
    if args.cumulant:
        script_name = os.path.join('cumulant_pipeline', 'single_cell_pipeline_cumulant.sh')
    elif args.impute:
        script_name = os.path.join('impute_pipeline', 'single_cell_pipeline_impute.sh')
    else:
        logging.info("Running standard single-cell pipeline")
        script_name = os.path.join('standard_pipeline', 'single_cell_pipeline.sh')

    # Construct full script path
    script_path = os.path.join(script_dir, 'utilities', 'single-cell', script_name)

    # Ensure the script is executable
    if not os.access(script_path, os.X_OK):
        logging.info(f"Adding executable permissions to {script_path}.")
        os.chmod(script_path, 0o755)

    # Execute the script
    run_shell_script(script_path)

def diffchrom(args):
    global config
    # Hi-C and epigenetic analysis based on the differentiate_chromosomes command
    hic_file = args.hic_file if args.hic_file else config.get('bam_directory')
    epigenetic_file = args.epigenetic_file if args.epigenetic_file else None
    output_dir = args.output_dir if args.output_dir else config.get('output_directory')

    logging.info(f"Running differentiate_chromosomes with Hi-C file: {hic_file}, Epigenetic file: {epigenetic_file}, Output directory: {output_dir}")
    # Placeholder for actual analysis code

def preprocess(args):
    global config
    logging.info("Preprocess BAM files for analysis")

    # Preprocess BAM files for analysis
    input_dir = args.input_dir
    output_dir = args.output_dir
    ref_genome = args.ref_genome

    # Update config with preprocess parameters
    updates = {
        'bam_directory': input_dir,
        'output_directory': output_dir,
        'reference_genome': ref_genome,
        'genome_fa': args.genome_fa,
        'chrom_file': args.chrom_sizes
    }
    config.update(updates)
    save_config(config)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set appropriate scripts based on the flags provided
    if args.nomehic:
        preprocess_script = 'nomehic_preprocess.sh'
    elif args.m3C:
        preprocess_script = 'm3C_preprocess.sh'
    elif args.test:
        preprocess_script = 'run_process_test.sh'
    else:
        logging.error("Either the --nomehic flag or the --m3C flag is required.")
        sys.exit(1)

    # Construct the full path to the preprocess script
    preprocess_script_path = os.path.join(script_dir, 'utilities', 'single-cell', 'preprocess', preprocess_script)

    # Ensure the script is executable
    if not os.access(preprocess_script_path, os.X_OK):
        logging.info(f"Adding executable permissions to {preprocess_script_path}.")
        os.chmod(preprocess_script_path, 0o755)

    # Execute the script
    run_shell_script(preprocess_script_path)

def run_shell_script(script_path):
    logging.info(f"Script path resolved to: {script_path}")

    # Ensure the script is executable
    if not os.access(script_path, os.X_OK):
        logging.info(f"Adding executable permissions to {script_path}.")
        os.chmod(script_path, 0o755)

    # Set environment variable for config path
    env = os.environ.copy()
    env['HYPERMATRIX_CONFIG'] = config_path

    # Execute the script
    command = [script_path]
    logging.info(f"Executing command: {' '.join(command)}")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
        logging.info(f"Command executed successfully. Output:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with return code {e.returncode}. Error message:\n{e.stderr}")
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()
