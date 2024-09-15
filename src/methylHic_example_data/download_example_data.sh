#/bin/bash

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
# Assuming a specific file from the extraction exists to verify that it's already been extracted
if [ ! -f "../../projects/methyHic/GSM3359999_JL_457_1.CGATGT.mm9.calmd.cpg.filtered.sort.CG.strand.6plus2.bed.gz" ]; then
    echo "Extracting GSE119171_RAW.tar..."
    tar -xvf ../../projects/methyHic/GSE119171_RAW.tar -C ../../projects/methyHic/
else
    echo "Files already extracted. Skipping extraction."
fi


# Use the manually created SRR_Acc_List.txt for downloading SRA data
while read SRR; do
    if [ ! -f "../projects/methyHic/${SRR}.sra" ]; then
        echo "Downloading ${SRR}..."
        prefetch -O ../projects/methyHic/ ${SRR}
    else
        echo "${SRR}.sra already exists. Skipping download."
    fi
done < SRR_Acc_List.txt

# Optionally, convert the downloaded SRA files to FASTQ
# Uncomment the following lines if you need FASTQ conversion
# echo "Converting SRA files to FASTQ..."
# for sra_file in ../projects/methyHic/*.sra; do
#     fastq-dump --split-files --gzip $sra_file -O ../projects/methyHic/
# done

echo "All downloads completed."

