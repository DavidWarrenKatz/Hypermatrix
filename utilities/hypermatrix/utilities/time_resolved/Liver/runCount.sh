#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=2-01:50:00
#SBATCH --ntasks=8
#SBATCH --mem=80000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetCount.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetCount.out

cd /home/dwk681/workspace/CRA004660/Liver

<<comment
/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403690_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Liver/CRR403690_Liver-Iso-Y \
                   --sample=CRR403690 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403691_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Liver/CRR403691_Liver-Iso-O \
                   --sample=CRR403691 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403692_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Liver/CRR403692_Liver-Het-Y \
                   --sample=CRR403692 \
                   --localcores=8 \
                   --localmem=64

comment

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR403693_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Liver/CRR403693_Liver-Het-O \
                   --sample=CRR403693 \
                   --localcores=8 \
                   --localmem=64

