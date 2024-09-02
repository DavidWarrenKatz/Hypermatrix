#!/usr/bin/env python3
# file: hypermatrix/main.py

import pkg_resources
import argparse
import subprocess
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))
config_path = os.path.join(script_dir)
sys.path.append(config_path)

# Import config from config_and_print
from config_and_print import config

VERSION = """
Hypermatrix version 0.1 - A tool for integrating
multi-omics data and epigenetic analysis using
advanced tensor techniques.
"""

def main():
    print("[DEBUG]: entering main function")
    parser = argparse.ArgumentParser(description='Hypermatrix command-line tool')
    parser.add_argument('-v','--version', action='version', version=VERSION)

    subparsers = parser.add_subparsers(dest='command')

    # Create the ABcluster subcommand
    abcluster_parser = subparsers.add_parser('ABcluster', help='Run the ABcluster analysis')
    abcluster_parser.add_argument('--methylation_file', type=str, help='Path to the single-cell CpG methylation file')
    abcluster_parser.add_argument('--conformation_file', type=str, help='Path to the chromosome conformation file')
    abcluster_parser.add_argument('--output_dir', type=str, help='Directory where the output files will be saved')
    abcluster_parser.add_argument('-cumulant', action='store_true', help='Run cumulant shell script')
    abcluster_parser.add_argument('-impute', action='store_true', help='Run impute shell script')

    # Create the differentiate_chromosomes subcommand
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
    else:
        print(f"Unknown command: {args.command}")
        parser.print_help()

def abcluster(args):
    methylation_file = args.methylation_file if args.methylation_file else config['methy_directory']
    conformation_file = args.conformation_file if args.conformation_file else config['bam_directory']
    output_dir = args.output_dir if args.output_dir else config['output_directory']

    # edit the following to point to the correct scripts
    # this update will ensure that the script is called from cite packages 
    # all packages will be called relative to room code directory 

    cumulant_script = r"single_cell_pipleline_cumulant.sh"
    if args.cumulant:
        run_shell_script(cumulant_script, methylation_file, conformation_file, output_dir)
    elif args.impute:
        run_shell_script('impute_script.sh', methylation_file, conformation_file, output_dir)
    else:
        print("No action specified. Use -cumulant or -impute to run a specific analysis.")

def diffchrom(args):
    hic_file = args.hic_file if args.hic_file else config['hic_directory']
    epigenetic_file = args.epigenetic_file if args.epigenetic_file else None
    output_dir = args.output_dir if args.output_dir else config['output_directory']

    # Placeholder for analysis code
    print(f"Running differentiate_chromosomes with Hi-C file: {hic_file}, Epigenetic file: {epigenetic_file}, Output directory: {output_dir}")
    # Add your analysis code here


def run_shell_script(script_name, methylation_file, conformation_file, output_dir):
    # Resolve the full path of the script using pkg_resources
    full_script_path = pkg_resources.resource_filename('hypermatrix', os.path.join('utilities', 'single-cell', script_name))
    print(f"[INFO]: Script path resolved to: {full_script_path}")

    # Ensure the script is executable
    if not os.access(full_script_path, os.X_OK):
        print(f"[INFO]: Adding executable permissions to {full_script_path}.")
        os.chmod(full_script_path, 0o755)

    # Display the command to be executed
    command = f"{full_script_path} {methylation_file} {conformation_file} {output_dir}"
    print(f"[INFO]: Executing command: {command}")

    # Run the command and capture output
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print(f"[INFO]: Command executed successfully. Output:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR]: Command failed with return code {e.returncode}. Error message:\n{e.stderr}")
        pass

if __name__ == "__main__":
    main()

