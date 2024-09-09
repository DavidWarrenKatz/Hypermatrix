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

# Import config file
from config import config

VERSION = """
Hypermatrix version 0.1 - A tool for integrating
multi-omics data and epigenetic analysis using
advanced tensor techniques.
"""

def main():
    #print("[DEBUG-1]: entering main function")
    parser = argparse.ArgumentParser(description='Hypermatrix command-line tool')
    parser.add_argument('-v','--version', action='version', version=VERSION)

    subparsers = parser.add_subparsers(dest='command')

    # Create the ABcluster subcommand
    # Perhaps add second shorter options later
    abcluster_parser = subparsers.add_parser('ABcluster', help='Run the ABcluster analysis')
    abcluster_parser.add_argument('--methy', type=str, help='Path to the single-cell CpG methylation directory')
    abcluster_parser.add_argument('--hic', type=str, help='Path to the single-cell chromosome conformation directory')
    abcluster_parser.add_argument('--output_dir', type=str, help='Directory where the output files will be saved')
    abcluster_parser.add_argument('-cumulant', action='store_true', help='Run cumulant shell script')
    abcluster_parser.add_argument('-impute', action='store_true', help='Run impute shell script')
    abcluster_parser.add_argument('-genome_id', type=str, help='Reference Genome')
    abcluster_parser.add_argument('-res', type=str, help='Resolution, Size of genomic bins')

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
        print(f"Unknown command: {args.command}. Try the ABcluster command")
        parser.print_help()

def update_config_file(config_file_path, new_methy_directory):
    """
    Function to update the methy_directory value in the config.py file.
    Args:
    - config_file_path (str): Path to the config.py file.
    - new_methy_directory (str): The new path for methy_directory.
    """
    with open(config_file_path, 'r') as file:
        lines = file.readlines()

    with open(config_file_path, 'w') as file:
        for line in lines:
            if line.startswith('methy_directory ='):
                file.write(f"methy_directory = '{new_methy_directory}'\n")
            else:
                file.write(line)

    print(f"[INFO]: Updated methy_directory in {config_file_path} to {new_methy_directory}")

def abcluster(args):
    # Overwrite config value if --methy is provided
    methy_directory = args.methy if args.methy else config['methy_directory']
    conformation_file = args.conformation_file if args.conformation_file else config['bam_directory']
    output_dir = args.output_dir if args.output_dir else config['output_directory']

    print("[DEBUG-12]: Resolving path to first script")
    cumulant_script = "/utilities/single-cell/single_cell_pipleline_cumulant.sh"
    impute_script = "path/to/script.sh"  # update later

    if args.cumulant:
        run_shell_script(cumulant_script, methylation_file, conformation_file, output_dir)
    elif args.impute:
        run_shell_script(impute_script, methylation_file, conformation_file, output_dir)
    else:
        print("No action specified. Use -cumulant or -impute to run a specific analysis.")


def abcluster(args):
    # Overwrite the config.py file if --methy is --hic or other flags are specified
    if args.methy:
        update_config_file(config_file_path, args.methy)
    
    standard_script = "./utilities/single-cell/single_cell_pipleline.sh" 
    cumulant_script = "./utilities/single-cell/single_cell_pipleline_cumulant.sh"
    impute_script = "./utilities/single-cell/single_cell_pipleline_impute.sh"  

    if args.cumulant:
        run_shell_script(cumulant_script, methylation_file, conformation_file, output_dir)
    elif args.impute:
        run_shell_script(impute_script, methylation_file, conformation_file, output_dir)
    else:
        run_shell_script(impute_script, methylation_file, conformation_file, output_dir)

def diffchrom(args):
    hic_file = args.hic_file if args.hic_file else config['hic_directory']
    epigenetic_file = args.epigenetic_file if args.epigenetic_file else None
    output_dir = args.output_dir if args.output_dir else config['output_directory']

    # Placeholder for analysis code
    print(f"Running differentiate_chromosomes with Hi-C file: {hic_file}, Epigenetic file: {epigenetic_file}, Output directory: {output_dir}")
    # Add your analysis code here

def next_fnction():
    pass

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
    # fails to execute bash command, clams file not found 
    # redudant
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print(f"[INFO]: Command executed successfully. Output:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR]: Command failed with return code {e.returncode}. Error message:\n{e.stderr}")
        pass

if __name__ == "__main__":
    main()

