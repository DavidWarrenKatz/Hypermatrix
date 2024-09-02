#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=4-05:50:00
#SBATCH --ntasks=8
#SBATCH --mem=80000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Brain/logs/%JgetCount.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Brain/logs/%JgetCount.out

cd /home/dwk681/workspace/CRA004660/Brain

<<comment
/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403686_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Brain/CRR403686_Brain-Iso-Y \
                   --sample=CRR403686 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403687_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Brain/CRR403687_Brain-Iso-O  \
                   --sample=CRR403687 \
                   --localcores=8 \
                   --localmem=64

comment

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403688_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Brain/CRR403688_Brain-Het-Y  \
                   --sample=CRR403688 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403689_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Brain/CRR403689_Brain-Het-O \
                   --sample=CRR403689 \
                   --localcores=8 \
                   --localmem=64













