##################################################################################
# generate the interval bed file used in the methylation processing for scNome-Seq
##################################################################################

import os
from config_and_print import chrom_file, resolutions, output_directory, software_directory

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
            size = int(size)
            for start in range(0, size, resolution):
                end = min(start + resolution, size)
                f_out.write(f"{chrom}\t{start}\t{end}\n")

# Define the output BED file path
output_bed_file = os.path.join(output_directory, f'hg19.common_chr.{resolution_label}_interval.autosome.bed')

# Check if the file already exists
if os.path.exists(output_bed_file):
    print(f"File '{output_bed_file}' already exists. Skipping the generation step.")
else:
    print(f"Generating interval BED file: '{output_bed_file}'")
    generate_interval_bed(chrom_file, resolution, output_bed_file)
    print(f"Interval BED file generated successfully: '{output_bed_file}'")

























