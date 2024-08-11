#!/bin/bash
# Description: The following script pulls bulk hi-c matrices for further processing downstream
# Planned upgrades 
# To-do: setup logging directory, capture script name, date and runtime 
# To-do: send all output to the logs file and also log to stdout by piping to tee 
# To-do: set up a small test mode by using only part of the bigwig files e.g. chr1&chr2, add -test param {test  mode}

# Exit immediately if a command exits with a non-zero status
set -e
# Logging function
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $@"
}

# Function to create directories
create_directory() {
    local dir_path=$1
    if [ ! -d "$dir_path" ]; then
        mkdir -p "$dir_path"
        log "Created directory: $dir_path"
    else
        log "Directory already exists: $dir_path"
    fi
}


# Function to download files if they don't exist
download_file_if_missing() {
    local file_path=$1
    local url=$2
    if [ ! -f "$file_path" ]; then
        wget "$url" -O "$file_path"
        log "Downloaded: $file_path"
    else
        log "File already exists: $file_path"
    fi
}

# Run the Python script and source the output to import the variables
log "Running config_and_print.py to load variables"
eval $(python3 config_and_print.py)

# Set up directories
create_directory "$data_path"
create_directory "${data_path}hicFiles/short_score_textform/"
create_directory "${data_path}Workspaces/individual/"

# Convert lists to strings for Python script
resolutions_list=$(printf "'%s', " "${resolutions[@]}" | sed 's/, $//')
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}" | sed 's/, $//')
data_types_list=$(printf "'%s', " "${data_types[@]}" | sed 's/, $//')

# Download necessary files
download_file_if_missing "../../bin/softwarefiles/dark_regions_hg19.bigWig" "$dark_regions_hg19_url"
download_file_if_missing "../../bin/softwarefiles/H3K4me3.bigwig" "$H3K4me3_GM12878_hg19_url"

# Execute the Python script to create dark regions file
log "Executing retrieve_dark_bins.py"
python retrieve_dark_bins.py

# Execute Python scripts for processing Hi-C data
log "Executing extract_hic_data.py"
python extract_hic_data.py "$data_path" "$hic_url" "$resolutions_list" "$chromosomes_list" "$data_types_list"

log "Executing process_hic_files_cumulant.py"
python process_hic_files_cumulant.py "$data_path" "$resolutions_list" "$chromosomes_list" "$data_types_list"

# Check and unzip tensorlab files if needed
zip_file="../../bin/softwarefiles/tensorlab/tensorlab_2016-03-28.zip"
unzip_dir="../../bin/softwarefiles/tensorlab/"

log "Checking for tensorlab files"
if [ -f "$zip_file" ]; then
    if all_files_exist; then
        log "All files from the zip are already present in the directory. Skipping unzip."
    else
        log "Unzipping tensorlab files..."
        unzip "$zip_file" -d "$unzip_dir"
    fi
else
    log "Tensorlab zip file does not exist. Exiting."
    exit 1
fi

# Load the necessary modules
if ! module list 2>&1 | grep -q "matlab"; then
    log "Loading MATLAB module"
    module load matlab/r2022b
else
    log "MATLAB module already loaded"
fi

# Execute the MATLAB script
log "Executing MATLAB script getStructuredData.m"
matlab -nodisplay -r "run('getStructuredData.m'); exit;"

# Run additional scripts
log "Setting permissions and executing scripts"
chmod +x check_and_run.sh run_java_commands.sh
./check_and_run.sh

log "Script completed successfully"


