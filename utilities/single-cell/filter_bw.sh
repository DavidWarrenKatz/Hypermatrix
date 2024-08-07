#!/bin/bash

# Load configuration using a Python script
declare -A config
while IFS='=' read -r key value; do
    # Remove leading and trailing quotes from the value
    value=$(echo "$value" | sed -e "s/^'//" -e "s/'$//")
    config[$key]=$value
done < <(python3 config_and_print.py)

# Load prefixes from filtered_bam_list.txt
prefixes=()
while IFS= read -r line; do
    prefixes+=("$line")
done < "${config[filtered_list]}"

# Print loaded configuration for debugging
echo "Loaded configuration:"
for key in "${!config[@]}"; do
    echo "$key=${config[$key]}"
done

# Configuration
java_path="/home/dwk681/workspace/softwareFiles/java/jdk1.8.0_281/bin/java"
java_opts="-Xmx18G"
java_classpath="${config[software_directory]}/dnaaseUtils-0.14-jar-with-dependencies.jar:${config[software_directory]}/java-genomics-io.jar:${config[software_directory]}/igv.jar"
main_class="main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed"
bed_file="${config[software_directory]}/b37.common_chr.1Mb_interval.autosome.bed"
output_bed_base="${config[output_directory]}/b37.autosome.1Mb_interval.add_value.methy.bed.gz"

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
    cov_file="${config[methy_directory]}/${prefix}.cov.b37.bw"
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
    printf "%s\n" "${not_found_prefixes[@]}" | grep -v '^[[:space:]]*$' | tee "${config[output_directory]}/missing_prefixes.txt"
else
    echo "All prefixes were found."
    rm -f "${config[output_directory]}/missing_prefixes.txt"  # Clean up the missing prefixes file if it's empty
fi

# Clean up
rm -f java_error.log methy_summary.cmd.txt





