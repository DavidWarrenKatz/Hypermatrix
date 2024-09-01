# hypermatrix/main.py

import os
import sys

# Determine the script's directory and add the path relative to it
script_dir = os.path.dirname(os.path.realpath(__file__))
config_path = os.path.join(script_dir)
sys.path.append(config_path)


import argparse
import subprocess

# Import config from config_and_print
from config_and_print import config

def main():
    parser = argparse.ArgumentParser(description='Hypermatrix command-line tool')
    subparsers = parser.add_subparsers(dest='command')

    # Create the ABcluster subcommand
    abcluster_parser = subparsers.add_parser('ABcluster', help='Run the ABcluster analysis')
    abcluster_parser.add_argument('--methylation_file', type=str, help='Path to the single-cell CpG methylation file')
    abcluster_parser.add_argument('--conformation_file', type=str, help='Path to the chromosome conformation file')
    abcluster_parser.add_argument('--output_dir', type=str, help='Directory where the output files will be saved')
    abcluster_parser.add_argument('-cumulant', action='store_true', help='Run cumulant shell script')
    abcluster_parser.add_argument('-impute', action='store_true', help='Run impute shell script')

    args = parser.parse_args()

    # Dispatch the command
    if args.command == 'ABcluster':
        abcluster(args)
    else:
        print(f"Unknown command: {args.command}")

def abcluster(args):
    # Override with provided arguments
    methylation_file = args.methylation_file if args.methylation_file else config['methy_directory']
    conformation_file = args.conformation_file if args.conformation_file else config['bam_directory']
    output_dir = args.output_dir if args.output_dir else config['output_directory']

    if args.cumulant:
        run_shell_script('cumulant_script.sh', methylation_file, conformation_file, output_dir)
    elif args.impute:
        run_shell_script('impute_script.sh', methylation_file, conformation_file, output_dir)
    else:
        print("No action specified. Use -cumulant or -impute to run a specific analysis.")

def run_shell_script(script_name, methylation_file, conformation_file, output_dir):
    # Replace placeholders in shell script with actual file paths
    command = f"./{script_name} {methylation_file} {conformation_file} {output_dir}"
    subprocess.run(command, shell=True)

if __name__ == "__main__":
    main()

