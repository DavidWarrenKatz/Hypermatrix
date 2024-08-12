from config_and_print import software_directory, resolutions, output_directory, filtered_list, normalization
import os
import subprocess
import time

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
juicer_tools_path = os.path.join(software_directory, "juicer_tools_1.22.01.jar")
prefix_file = filtered_list
normalization_method = normalization  # Set to "KR" for Knight-Ruiz or "ICE" for Iterative Correction

# List of autosomal chromosomes
chromosomes = [str(i) for i in range(1, 23)]  # Chromosomes 1 to 22

# Phase 1: Create all .hic files in a dictionary
def create_all_hic_files():
    hic_files = {}
    for chromosome in chromosomes:
        for dataset in datasets:
            hic_file = create_hic_file(chromosome, dataset)
            hic_files[(chromosome, dataset)] = hic_file
    return hic_files

def convert_to_hic_short_format(input_file, chrom, output_file, resolution):
    # Check if the formatted file already exists
    if os.path.exists(output_file):
        print(f"Formatted file already exists: {output_file}. Skipping conversion.")
        return  # Skip the conversion if the file exists
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) != 3:
                print(f"Skipping invalid line in {input_file}: {line.strip()}")
                continue
            
            # Extract the values from the short form
            chr1 = chr2 = chrom  # both positions belong to the same chromosome
            pos1 = int(parts[0]) * resolution
            pos2 = int(parts[1]) * resolution
            score = parts[2]

            # Construct the HiC Short format line
            str1 = "0"
            str2 = "0"
            frag1 = "0"
            frag2 = "1"

            # Writing the line in HiC Short format with score
            outfile.write(f"{str1} {chr1} {pos1} {frag1} {str2} {chr2} {pos2} {frag2} {score}\n")

    print(f"Conversion completed: {output_file}")

# Function to create .hic file from short format with the selected normalization method
def create_hic_file(chromosome, dataset):
    short_form_file = os.path.join(data_path, f"chr{chromosome}", f"{dataset}_chr{chromosome}.txt")
    formatted_file = os.path.join(data_path, f"chr{chromosome}", f"{dataset}_chr{chromosome}_juicer_formatted.txt")
    hic_file = os.path.join(data_path, f"chr{chromosome}", f"{dataset}_chr{chromosome}.hic")

    # Convert to HiC Short format with score using the provided resolution
    convert_to_hic_short_format(short_form_file, chromosome, formatted_file, resolution)

    # Proceed with .hic file creation
    if not is_valid_hic_file(hic_file):
        java_command = [
            "java", "-Xmx100G", "-jar", juicer_tools_path, "pre",
            "-r", str(resolution), formatted_file, hic_file, "hg19"
        ]
        print(f"Running: {' '.join(java_command)}")
        
        try:
            subprocess.run(java_command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during HIC file creation: {e}")
            return None  # Return None to indicate failure
    
    # Retry mechanism
    max_retries = 10
    for i in range(max_retries):
        if is_valid_hic_file(hic_file) and validate_hic_file(hic_file, chromosome):
            print(f"HIC file successfully created: {hic_file}")
            return hic_file
        else:
            print(f"HIC file not created yet or is invalid, retrying ({i + 1}/{max_retries})...")
            time.sleep(20)  # Wait 10 seconds before retrying
    
    print(f"Failed to create a valid HIC file after {max_retries} retries: {hic_file}")
    return None

# Function to check if the .hic file is valid based on file size and validation command
def is_valid_hic_file(hic_file, min_size=10):
    return os.path.isfile(hic_file) and os.path.getsize(hic_file) > min_size

# Function to validate .hic file using juicer_tools
def validate_hic_file(hic_file, chrom):
    validate_command = [
        "java", "-jar", juicer_tools_path, "dump", "observed", "NONE",
        hic_file, chrom, chrom, "BP", str(resolution), "/dev/null"
    ]
    try:
        subprocess.run(validate_command, check=True)
        return True
    except subprocess.CalledProcessError:
        return False

# Phase 2: Apply KR normalization to all created .hic files
def apply_kr_normalization(hic_files):
    for (chromosome, dataset), hic_file in hic_files.items():
        kr_normalize_and_convert(hic_file, chromosome, dataset)

# Function to apply KR normalization and convert back to short format
def kr_normalize_and_convert(hic_file, chromosome, dataset):
    if hic_file is None:
        print(f"Skipping KR normalization due to missing or invalid HIC file for {dataset} chr{chromosome}.")
        return

    kr_matrix_file = os.path.join(data_path, f"chr{chromosome}", f"{dataset}_chr{chromosome}_KR_matrix.txt")
    kr_short_file = os.path.join(data_path, f"chr{chromosome}", f"{dataset}_chr{chromosome}_KR_short.txt")

    java_command = [
        "java", "-jar", juicer_tools_path, "dump", "observed", "KR",
        hic_file, f"chr{chromosome}", f"chr{chromosome}", "BP", str(resolution), kr_matrix_file
    ]
    print(f"Running KR normalization: {' '.join(java_command)}")
    
    try:
        subprocess.run(java_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during KR normalization: {e}")
        return

    with open(kr_matrix_file, 'r') as f_in, open(kr_short_file, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) == 3:
                bin1, bin2, count = parts
                f_out.write(f"{bin1}\t{bin2}\t{count}\n")
    
    print(f"KR normalization completed for {dataset} chr{chromosome}. Short format saved to {kr_short_file}")

# Main execution flow
# Read prefixes from the file
with open(prefix_file, 'r') as f:
    datasets = f.read().splitlines()

# Phase 1: Create all .hic files
hic_files = create_all_hic_files()

# Phase 2: Apply KR normalization to all created .hic files
apply_kr_normalization(hic_files)

print("Processing completed.")

