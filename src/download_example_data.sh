#!/bin/bash

# Load SRA Toolkit module
module load sra-toolkit/2.10.9

# Create directory if it doesn't exist
mkdir -p ../projects/methyHic

# Check if the main tar file is already downloaded
if [ ! -f ../projects/methyHic/GSE119171_RAW.tar ]; then
    echo "Downloading GSE119171_RAW.tar..."
    wget -O ../projects/methyHic/GSE119171_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119171/suppl/GSE119171_RAW.tar"
else
    echo "GSE119171_RAW.tar already exists. Skipping download."
fi

# Check if the files from the tar archive are already extracted
# Assuming a specific file from the extraction exists to verify that it's already been extracted
if [ ! -f "../projects/methyHic/GSM3359999_JL_457_1.CGATGT.mm9.calmd.cpg.filtered.sort.CG.strand.6plus2.bed.gz" ]; then
    echo "Extracting GSE119171_RAW.tar..."
    tar -xvf ../projects/methyHic/GSE119171_RAW.tar -C ../projects/methyHic/
else
    echo "Files already extracted. Skipping extraction."
fi

#create a file SRR_Acc_List.txt
'''
SRR7770795
SRR7770796
SRR7770797
SRR7770798
SRR7770799
SRR7770800
SRR7770801
SRR7770802
SRR7770803
SRR7770804
SRR7770805
SRR7770806
SRR7770807
SRR7770808
SRR7770809
SRR7770810
SRR7770811
SRR7770812
SRR7770813
SRR7770814
SRR7770815
SRR7770816
SRR7770817
SRR7770818
SRR7770819
SRR7770820
SRR7770821
SRR7770822
SRR7770823
SRR7770824
SRR7770825
SRR7770826
SRR7770827
SRR7770828
SRR7770829
SRR7770830
SRR7770831
SRR7770832
SRR7770833
SRR7770834
SRR7770835
SRR7770836
SRR7770837
SRR7770838
SRR7770839
SRR7770840
SRR7770841
SRR7770842
SRR7770843
SRR7770844
SRR7770845
SRR7770846
SRR7770847
SRR7770848
SRR7770849
SRR7770850
SRR7770851
SRR7770852
SRR7770853
SRR7770854
SRR7770855
SRR7770856
SRR7770857
SRR7770858
SRR7770859
SRR7770860
SRR7770861
SRR7770862
SRR7770863
SRR7770864
SRR7770865
SRR7770866
SRR7770867
SRR7770868
SRR7770869
SRR7770870
SRR7770871
SRR7770872
SRR7770873
SRR7770874
SRR7770875
SRR7770876
SRR7770877
SRR7770878
SRR7770879
SRR7770880
'''



# Download PRJNA488313 data from SRA
# Fetch all the SRR numbers associated with PRJNA488313 using the SRA Toolkit

if [ ! -f ../projects/methyHic/PRJNA488313_sra_list.txt ]; then
    echo "Fetching SRR list for PRJNA488313..."
    esearch -db sra -query PRJNA488313 | efetch -format runinfo | cut -d ',' -f1 | grep SRR > ../projects/methyHic/PRJNA488313_sra_list.txt
else
    echo "SRR list already exists. Skipping fetching."
fi

# Loop through each SRR in the list and download if not already present
while read SRR; do
    if [ ! -f "../projects/methyHic/${SRR}.sra" ]; then
        echo "Downloading ${SRR}..."
        prefetch -O ../projects/methyHic/ ${SRR}
    else
        echo "${SRR}.sra already exists. Skipping download."
    fi
done < ../projects/methyHic/PRJNA488313_sra_list.txt

# Optionally, convert the downloaded SRA files to FASTQ
# Uncomment the following lines if you need FASTQ conversion
# echo "Converting SRA files to FASTQ..."
# for sra_file in ../projects/methyHic/*.sra; do
#     fastq-dump --split-files --gzip $sra_file -O ../projects/methyHic/
# done

echo "All downloads completed."

