#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=4-01:50:00
#SBATCH --ntasks=8
#SBATCH --mem=80000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Skin/logs/%JgetCount.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Skin/logs/%JgetCount.out

cd /home/dwk681/workspace/CRA004660/Skin

<<comment
/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308456_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308456_Skin-Iso-Y-1 \
                   --sample=CRR308456 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308457_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308457_Skin-Iso-Y-2 \
                   --sample=CRR308457 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308458_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308458_Skin-Iso-Y-3\
                   --sample=CRR308458 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308459_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308459_Skin-Iso-Y-4 \
                   --sample=CRR308459 \
                   --localcores=8 \
                   --localmem=64


/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308460_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308460_Skin-Iso-Y-5 \
                   --sample=CRR308460 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308461_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308461_Skin-Iso-Y-6\
                   --sample=CRR308461 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308462_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/Iso-Y/CRR308462_Skin-Iso-Y-7 \
                   --sample=CRR308462 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308463_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/Iso-Y/CRR308463_Skin-Iso-Y-8 \
                   --sample=CRR308463 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308464_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308464_Skin-Iso-O-1 \
                   --sample=CRR308464 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308465_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308465_Skin-Iso-O-2 \
                   --sample=CRR308465 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308466_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308466_Skin-Iso-O-3 \
                   --sample=CRR308466 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308467_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308467_Skin-Iso-O-4 \
                   --sample=CRR308467 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308468_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308468_Skin-Iso-O-5 \
                   --sample=CRR308468 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308469_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308469_Skin-Iso-O-6 \
                   --sample=CRR308469 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308470_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308470_Skin-Iso-O-7 \
                   --sample=CRR308470 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308471_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308471_Skin-Iso-O-8 \
                   --sample=CRR308471 \
                   --localcores=8 \
                   --localmem=64


/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308472_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308472_Skin-Het-Y-1 \
                   --sample=CRR308472 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308473_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308473_Skin-Het-Y-2 \
                   --sample=CRR308473 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308474_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308474_Skin-Het-Y-3 \
                   --sample=CRR308474 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308475_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308475_Skin-Het-Y-4 \
                   --sample=CRR308475 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308476_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308476_Skin-Het-O-1 \
                   --sample=CRR308476 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308477_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308477_Skin-Het-O-2 \
                   --sample=CRR308477 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308478_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308478_Skin-Het-O-3 \
                   --sample=CRR308478 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308479_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308479_Skin-Het-O-4 \
                   --sample=CRR308479 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308480_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308480_Skin-Het-O-5 \
                   --sample=CRR308480 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308481_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308481_Skin-Het-O-6 \
                   --sample=CRR308481 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308482_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308482_Skin-Het-O-7 \
                   --sample=CRR308482 \
                   --localcores=8 \
                   --localmem=64

/home/dwk681/workspace/softwareFiles/cellranger-7.2.0/./cellranger count --id=CRR308483_cellranger_output \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/Skin/CRR308483_Skin-Het-O-8 \
                   --sample=CRR308483 \
                   --localcores=8 \
                   --localmem=64

comment

