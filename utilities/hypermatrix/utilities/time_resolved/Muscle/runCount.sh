#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=4-05:50:00
#SBATCH --ntasks=8
#SBATCH --mem=80000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Muscle/logs/%JgetCount.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Muscle/logs/%JgetCount.out

cd /home/dwk681/workspace/CRA004660/Muscle

<<comment
/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403694_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Muscle/CRR403694_Muscle-Iso-Y \
                   --sample=CRR403694 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403696_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Muscle/CRR403696_Muscle-Het-Y  \
                   --sample=CRR403696 \
                   --localcores=8 \
                   --localmem=64


/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR477130_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Muscle/CRR477130_Muscle-Het-Y  \
                   --sample=CRR477130 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403695_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Muscle/CRR403695_Muscle-Iso-O \
                   --sample=CRR403695 \
                   --localcores=8 \
                   --localmem=64
comment

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403697_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Muscle/CRR403697_Muscle-Het-O \
                   --sample=CRR403697 \
                   --localcores=8 \
                   --localmem=64










