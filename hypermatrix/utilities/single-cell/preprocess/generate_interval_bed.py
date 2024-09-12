##################################################################################
# generate the interval bed file used in the methylation processing for scNome-Seq
##################################################################################

import sys
import os

# Add the directory where config.py is located to the Python path
config_dir = '../../../'
config_dir = os.path.abspath(config_dir)  # Get absolute path

config_file = os.path.join(config_dir, 'config.py')  # Full path to config.py

# Check if the directory and config.py exist
if os.path.isdir(config_dir) and os.path.isfile(config_file):
    sys.path.append(config_dir)
    print(f"Config directory added to sys.path: {config_dir}")
    print(f"Found config.py at: {config_file}")
else:
    raise FileNotFoundError(f"config.py not found in directory: {config_dir}")

from config import chrom_file, resolutions, output_directory, software_directory, reference_genome

# Ensure resolutions is treated as a tuple or list of strings
if isinstance(resolutions, str):
    resolutions = (resolutions,)

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]

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

resolution, resolution_label = parse_resolution(resolution_str)

def generate_interval_bed(chrom_size_file, resolution, output_bed_file):
    with open(chrom_size_file, 'r') as f:
        chromosome_sizes = [line.strip().split() for line in f.readlines()]

    with open(output_bed_file, 'w') as f_out:
        for chrom, size in chromosome_sizes:
            # Remove "chr" prefix
            chrom = chrom.replace('chr', '')
            size = int(size)
            for start in range(0, size, resolution):
                end = min(start + resolution, size)
                f_out.write(f"{chrom}\t{start}\t{end}\n")

# Define the output BED file path
output_bed_file = os.path.join(software_directory, f'{reference_genome}.common_chr.{resolution_label}_interval.autosome.bed')

# Check if the file already exists
if os.path.exists(output_bed_file):
    print(f"File '{output_bed_file}' already exists. Skipping the generation step.")
else:
    print(f"Generating interval BED file: '{output_bed_file}'")
    generate_interval_bed(chrom_file, resolution, output_bed_file)
    print(f"Interval BED file generated successfully: '{output_bed_file}'")

