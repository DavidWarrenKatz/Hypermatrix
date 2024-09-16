#!/bin/bash

# Load SRA Toolkit module
module load sra-toolkit/2.10.9

# Create directory if it doesn't exist
mkdir -p ../../projects/methyHic

# Check if the main tar file is already downloaded
if [ ! -f ../../projects/methyHic/GSE119171_RAW.tar ]; then
    echo "Downloading GSE119171_RAW.tar..."
    wget -O ../../projects/methyHic/GSE119171_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119171/suppl/GSE119171_RAW.tar"
else
    echo "GSE119171_RAW.tar already exists. Skipping download."
fi

# Check if the files from the tar archive are already extracted
if [ ! -f "../../projects/methyHic/GSM3359999_JL_457_1.CGATGT.mm9.calmd.cpg.filtered.sort.CG.strand.6plus2.bed.gz" ]; then
    echo "Extracting GSE119171_RAW.tar..."
    tar -xvf ../../projects/methyHic/GSE119171_RAW.tar -C ../../projects/methyHic/
else
    echo "Files already extracted. Skipping extraction."
fi

# Use the manually created SRR_Acc_List.txt for downloading SRA data
while read SRR; do
    # Check if the directory for the SRR ID exists and skip download if it does
    if [ ! -d "../../projects/methyHic/${SRR}" ]; then
        echo "Downloading ${SRR}..."
        prefetch -O ../../projects/methyHic/ ${SRR}
    else
        echo "${SRR}.sra already exists in ${SRR} directory. Skipping download."
    fi
done < SRR_Acc_List.txt

# Convert all SRA files in each SRR directory to FASTQ format
for sra_dir in ../../projects/methyHic/SRR*/; do
    # Count the number of .sra files in the directory
    file_count=$(ls -1 "$sra_dir"/*.sra 2>/dev/null | wc -l)

    # Print the number of .sra files
    echo "Found $file_count SRA file(s) in $sra_dir"

    # Proceed with conversion if .sra files are found
    if [ $file_count -gt 0 ]; then
        for sra_file in "$sra_dir"/*.sra; do
            if [ -f "$sra_file" ]; then
                # Get the base name of the SRA file (without extension)
                base_name=$(basename "$sra_file" .sra)

                # Check if FASTQ files already exist
                if [ -f "../../projects/methyHic/${base_name}_1.fastq.gz" ] && [ -f "../../projects/methyHic/${base_name}_2.fastq.gz" ]; then
                    echo "FASTQ files for $sra_file already exist. Skipping conversion."
                else
                    echo "Converting $sra_file to FASTQ format..."
                    fastq-dump --split-files --gzip "$sra_file" -O ../../projects/methyHic/
                fi
            fi
        done
    else
        echo "No SRA files found in $sra_dir."
    fi
done

echo "All downloads and conversions completed."



