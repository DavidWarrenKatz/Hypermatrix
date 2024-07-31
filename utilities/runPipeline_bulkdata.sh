#!/bin/bash
# Description: The following script pulls bulk hi-c matrices for further processing downstream

# Variables
# Run the Python script and source the output to import the variables
eval $(python3 config_and_print.py)

# Set the path for data storage
mkdir -p $data_path # mkdir if it doesn't exist

# Convert lists to strings for Python script
resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }

short_form_directory="${data_path}hicFiles/short_score_textform/"
mkdir -p $short_form_directory # mkdir if it doesn't exist
            
individual_workspace_directory="${data_path}Workspaces/individual/"
mkdir -p $individual_workspace_directory # mkdir if it doesn't exist

# Download the dark regions file if it doesn't already exist
dark_regions_file="../projects/softwarefiles/ENCFF000EHJ_hg19_wgEncodeCrgMapabilityAlign36mer.bigWig"
if [ ! -f "$dark_regions_file" ]; then
    wget https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig -O "$dark_regions_file"
fi

# Execute the Python script to create dark regions file
python retrieve_dark_bins.py

# Execute the Python script
python extract_hic_data.py "$data_path" "$hic_url" "$resolutions_list" "$chromosomes_list" "$data_types_list"

python process_hic_files.py "$data_path" "$resolutions_list" "$chromosomes_list" "$data_types_list"

#Check to make sure tensorlab is available
# Define the paths
zip_file="../projects/softwarefiles/tensorlab/tensorlab_2016-03-28.zip"
unzip_dir="../projects/softwarefiles/tensorlab/"

# Check if the zip file exists
if [ ! -f "$zip_file" ]; then
  echo "Zip file does not exist."
  exit 1
fi

# Get the list of files in the zip archive
zip_files=$(unzip -l "$zip_file" | awk 'NR>3 {print $4}' | sed '$d')

# Function to check if all files from the zip are already in the directory
all_files_exist() {
  for file in $zip_files; do
    if [ ! -f "$unzip_dir/$file" ]; then
      return 1
    fi
  done
  return 0
}

# Check if all files exist
if all_files_exist; then
  echo "All files from the zip are already present in the directory. Skipping unzip."
else
  echo "Unzipping the file..."
  unzip "$zip_file" -d "$unzip_dir"
  echo "Unzipping completed."
fi

# Load the necessary modules
module load matlab/r2022b

# Execute the MATLAB script
matlab -nodisplay -r "run('getStructuredData.m'); exit;"

chmod +x check_and_run.sh
chmod +x run_java_commands.sh
./check_and_run.sh


