#!/bin/bash                                                                                                                                                                                                                                     
eval "$(python3 config_and_print.py)"

# Iterate over files matching the pattern "sc*.hic.txt" in the output directory
for in_file in $output_directory/sc*.hic.txt; do
  # Check if the file exists
  if [[ -f "$in_file" ]]; then
    # Prepare the output filename by replacing ".hic.txt" with ".hic_matrix.txt.gz"
    out_file="${in_file/.hic.txt/.hic_matrix.txt.gz}"

    # Construct a sanitized job name to avoid issues with special characters
    job_name="process_${in_file//[^a-zA-Z0-9]/_}"

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

