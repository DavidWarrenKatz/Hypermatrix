#filter_bw.sh
#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Import the parameters from config.py (relative to the script's directory)
eval "$(python3 "$SCRIPT_DIR/../../../export_config.py")"

# Extract resolution and label from the resolutions string
IFS=':' read -r resolution resolution_label <<< "$resolutions"

# Load prefixes from filtered_bam_list.txt
prefixes=()
while IFS= read -r line; do
    prefixes+=("$line")
done < "$filtered_list"

# Print loaded configuration for debugging
echo "Loaded configuration:"
echo "bam_directory=$bam_directory"
echo "methy_directory=$methy_directory"
echo "software_directory=$software_directory"
echo "output_directory=$output_directory"
echo "resolutions=$resolutions"
echo "filtered_list=$filtered_list"

# Dynamically find the most appropriate version of Java
java_path=$(command -v java)
if [[ -z "$java_path" ]]; then
    echo "Java not found. Please install Java and try again."
    exit 1
fi

java_opts="-Xmx18G"
java_classpath="$software_directory/dnaaseUtils-0.14-jar-with-dependencies.jar:$software_directory/java-genomics-io.jar:$software_directory/igv.jar"
main_class="main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed"
bed_file="$software_directory/${reference_genome}.common_chr.${resolution_label}_interval.autosome.bed"
output_bed_base="$output_directory/${reference_genome}.autosome.${resolution_label}_interval.add_value.methy.bed.gz"

# Check if the desired output file already exists
if [[ -f "$output_bed_base" ]]; then
    echo "Output file $output_bed_base already exists. Skipping computation."
    exit 0
fi

# Open the output file for writing command strings
output_cmd_file="methy_summary.cmd.txt"
: > "$output_cmd_file"  # Create or clear the file if it exists

# Track prefixes of files that were not found
not_found_prefixes=()

# Process each prefix to create the command strings
for prefix in "${prefixes[@]}"; do
    cov_file="$methy_directory/${prefix}.cov.b37.bw"
    methy_count_file="${cov_file/cov/methy_count}"
    if [[ -f "$cov_file" ]]; then
        echo " -bigWig $methy_count_file -useMean0 0 -regionMode 0 -bigWig $cov_file -useMean0 0 -regionMode 0" >> "$output_cmd_file"
    else
        not_found_prefixes+=("$prefix")
    fi
done

# Read the command strings from the file and join them into a single line
cmd=$(tr '\n' ' ' < "$output_cmd_file")

# Construct the full Java command
full_cmd="$java_path $java_opts -cp \"$java_classpath\" $main_class $bed_file $output_bed_base $cmd"

# Print the command to be executed
echo "Executing command: $full_cmd"

# Execute the Java command directly and capture error output
eval "$full_cmd 2>java_error.log"

# Check if the command was successful
if [ $? -ne 0 ]; then
    echo "Failed to execute Java command"
    cat java_error.log
    rm -f $output_bed_base
    exit 1
fi

echo "Computation completed successfully"

# Report and save prefixes of files that were not found, removing empty lines
if [ ${#not_found_prefixes[@]} -ne 0 ]; then
    echo "The following prefixes were not found:"
    printf "%s\n" "${not_found_prefixes[@]}" | grep -v '^[[:space:]]*$' | tee "$output_directory/missing_prefixes.txt"
else
    echo "All prefixes were found."
    rm -f "$output_directory/missing_prefixes.txt"  # Clean up the missing prefixes file if it's empty
fi

# Clean up
rm -f java_error.log methy_summary.cmd.txt

