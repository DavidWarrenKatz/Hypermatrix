from config_and_print import software_directory, resolutions, output_directory, filtered_list, normalization
import os
import subprocess

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

data_path = f"{output_directory}/hic_{resolution_label}_raw_dir"  # Base data path
juicer_tools_path = "../../bin/softwarefiles/juicer_tools_1.22.01.jar"
refix_file = filtered_list
normalization_method = normalization  # Set to "KR" for Knight-Ruiz or "ICE" for Iterative Correction

# List of autosomal chromosomes
chromosomes = [str(i) for i in range(1, 23)]  # Chromosomes 1 to 22

# Function to create .hic file from short format with the selected normalization method
def create_hic_file(chromosome, dataset):
    short_form_file = os.path.join(data_path, chromosome, f"{dataset}_chr{chromosome}.txt")
    hic_file = os.path.join(data_path, chromosome, f"{dataset}_chr{chromosome}.hic")

    # Check if the .hic file already exists
    if not os.path.isfile(hic_file):
        java_command = [
            "java", "-Xmx100G", "-jar", juicer_tools_path, "pre",
            "-r", resolution, short_form_file, hic_file, "hg19"
        ]
        print(f"Running: {' '.join(java_command)}")
        subprocess.run(java_command, check=True)
    else:
        print(f"HIC file already exists: {hic_file}")

# Read prefixes from the file
with open(prefix_file, 'r') as f:
    datasets = f.read().splitlines()

# Loop through autosomal chromosomes and datasets to process
for chromosome in chromosomes:
    for dataset in datasets:
        create_hic_file(chromosome, dataset)

print("Processing completed.")

