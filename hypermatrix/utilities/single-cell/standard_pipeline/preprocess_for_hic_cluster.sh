#preprocess_for_hic_cluster.sh

######################################################################################                   
#This shell script takes the hic.txt files and puts them into the form scHICluster
#expects.
#Consider rewriting this and the next part so that we do not use scHIcluster initially
#Author: David Katz (davidkatz02@gmail.com)                                                              #####################################################################################

#!/bin/bash                   

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Import the parameters from config.py (relative to the script's directory)
eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

# Iterate over files matching the pattern "sc*.hic.txt" in the output directory
for in_file in $output_directory/sc*.hic.txt; do
  # Skip files containing ".good_reads.hic.txt"
  if [[ "$in_file" == *".good_reads.hic.txt" ]]; then
    continue
  fi

  # Check if the file exists
  if [[ -f "$in_file" ]]; then
    # Prepare the output filename by replacing ".hic.txt" with ".hic_matrix.txt.gz"
    out_file="${in_file/.hic.txt/.hic_matrix.txt.gz}"

    # Check if the output file already exists
    if [[ -f "$out_file" ]]; then
      echo "Output file $out_file already exists, skipping."
      continue
    fi

    # Print the job information
    echo "Processing $in_file to $out_file"

    # Execute the preprocessing steps
    cut -d' ' -f1-7 "$in_file" | \
    sed 's/ /\t/g' | \
    perl -ne 'chomp; @f=split "\t"; $f[1]="chr$f[1]"; $f[5]="chr$f[5]"; print join("\t", @f)."\n";' | \
    gzip -c > "$out_file"

    # Check if the command was successful
    if [ $? -eq 0 ]; then
      echo "Successfully processed $in_file to $out_file"
    else
      echo "Failed to process $in_file" >&2
    fi
  else
    echo "No files matching pattern found."
  fi
done


