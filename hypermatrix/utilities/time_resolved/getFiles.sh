#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=8
#SBATCH --mem=80000
#SBATCH --error=/home/dwk681/workspace/CRA004660/logs/%JgetFiles.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/logs/%JgetFiles.out

<<comment
mkdir Muscle
cd Muscle
#Muscle

mkdir CRR477130_Muscle-Het-Y
cd CRR477130_Muscle-Het-Y
#wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR477130/CRR477130_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR477130/CRR477130_r2.fq.gz 
cd ..

mkdir CRR403696_Muscle-Het-Y
cd CRR403696_Muscle-Het-Y
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403696/CRR403696_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403696/CRR403696_r2.fq.gz
cd ..

mkdir CRR403697_Muscle-Het-O
cd CRR403697_Muscle-Het-O
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403697/CRR403697_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403697/CRR403697_r2.fq.gz
cd ..

mkdir CRR403694_Muscle-Iso-Y
cd CRR403694_Muscle-Iso-Y
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403694/CRR403694_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403694/CRR403694_r2.fq.gz
cd ..

mkdir CRR403695_Muscle-Iso-O
cd CRR403695_Muscle-Iso-O
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403695/CRR403695_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403695/CRR403695_r2.fq.gz
cd ..



mkdir Brain
cd Brain

mkdir CRR403688_Brain-Het-Y
cd CRR403688_Brain-Het-Y
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403688/CRR403688_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403688/CRR403688_r2.fq.gz
cd ..

mkdir CRR403686_Brain-Iso-Y
cd CRR403686_Brain-Iso-Y
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403686/CRR403686_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403686/CRR403686_r2.fq.gz
cd ..

mkdir CRR403687_Brain-Iso-O
cd CRR403687_Brain-Iso-O
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403687/CRR403687_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403687/CRR403687_r2.fq.gz
cd ..

mkdir CRR403689_Brain-Het-O
cd CRR403689_Brain-Het-O
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403689/CRR403689_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403689/CRR403689_r2.fq.gz
cd ..
comment

#mkdir Skin
cd Skin

<<comment
mkdir CRR308456_Skin-Iso-Y-1 
cd CRR308456_Skin-Iso-Y-1
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308456/CRR308456_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308456/CRR308456_r2.fq.gz
cd ..

mkdir CRR308457_Skin-Iso-Y-2
cd CRR308457_Skin-Iso-Y-2
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308457/CRR308457_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308457/CRR308457_r2.fq.gz
cd ..

mkdir CRR308458_Skin-Iso-Y-3
cd CRR308458_Skin-Iso-Y-3
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308458/CRR308458_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308458/CRR308458_r2.fq.gz
cd ..

mkdir CRR308459_Skin-Iso-Y-4
cd CRR308459_Skin-Iso-Y-4
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308459/CRR308459_f1.fq.gz 
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308459/CRR308459_r2.fq.gz
cd ..

mkdir CRR308460_Skin-Iso-Y-5
cd CRR308460_Skin-Iso-Y-5
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308460/CRR308460_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308460/CRR308460_r2.fq.gz
cd ..

mkdir CRR308461_Skin-Iso-Y-6
cd CRR308461_Skin-Iso-Y-6
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308461/CRR308461_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308461/CRR308461_r2.fq.gz
cd ..

mkdir CRR308462_Skin-Iso-Y-7
cd CRR308462_Skin-Iso-Y-7
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308462/CRR308462_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308462/CRR308462_r2.fq.gz
cd ..

mkdir CRR308463_Skin-Iso-Y-8
cd CRR308463_Skin-Iso-Y-8
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308463/CRR308463_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308463/CRR308463_r2.fq.gz
cd ..


mkdir CRR308464_Skin-Iso-O-1
cd CRR308464_Skin-Iso-O-1 
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308464/CRR308464_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308464/CRR308464_r2.fq.gz
cd ..

mkdir CRR308465_Skin-Iso-O-2
cd CRR308465_Skin-Iso-O-2
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308465/CRR308465_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308465/CRR308465_r2.fq.gz
cd ..

mkdir CRR308466_Skin-Iso-O-3
cd CRR308466_Skin-Iso-O-3
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308466/CRR308466_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308466/CRR308466_r2.fq.gz
cd ..

mkdir CRR308467_Skin-Iso-O-4
cd CRR308467_Skin-Iso-O-4
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308467/CRR308467_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308467/CRR308467_r2.fq.gz
cd ..

mkdir CRR308468_Skin-Iso-O-5
cd CRR308468_Skin-Iso-O-5
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308468/CRR308468_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308468/CRR308468_r2.fq.gz
cd ..

mkdir CRR308469_Skin-Iso-O-6
cd CRR308469_Skin-Iso-O-6
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308469/CRR308469_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308469/CRR308469_r2.fq.gz
cd ..

mkdir CRR308470_Skin-Iso-O-7
cd CRR308470_Skin-Iso-O-7
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308470/CRR308470_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308470/CRR308470_r2.fq.gz
cd ..

mkdir CRR308471_Skin-Iso-O-8
cd CRR308471_Skin-Iso-O-8
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308471/CRR308471_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308471/CRR308471_r2.fq.gz
cd ..

mkdir CRR308472_Skin-Het-Y-1
cd CRR308472_Skin-Het-Y-1
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308472/CRR308472_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308472/CRR308472_r2.fq.gz
cd ..

mkdir CRR308473_Skin-Het-Y-2
cd CRR308473_Skin-Het-Y-2
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308473/CRR308473_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308473/CRR308473_r2.fq.gz
cd ..

comment

#mkdir CRR308474_Skin-Het-Y-3
cd CRR308474_Skin-Het-Y-3
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308474/CRR308474_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308474/CRR308474_r2.fq.gz
cd ..


mkdir CRR308475_Skin-Het-Y-4
cd CRR308475_Skin-Het-Y-4
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308475/CRR308475_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308475/CRR308475_r2.fq.gz
cd ..

mkdir CRR308476_Skin-Het-O-1
cd CRR308476_Skin-Het-O-1
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308476/CRR308476_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308476/CRR308476_r2.fq.gz
cd ..

mkdir CRR308477_Skin-Het-O-2
cd CRR308477_Skin-Het-O-2
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308477/CRR308477_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308477/CRR308477_r2.fq.gz
cd ..

mkdir CRR308478_Skin-Het-O-3
cd CRR308478_Skin-Het-O-3
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308478/CRR308478_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308478/CRR308478_r2.fq.gz
cd ..

mkdir CRR308479_Skin-Het-O-4
cd CRR308479_Skin-Het-O-4
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308479/CRR308479_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308479/CRR308479_r2.fq.gz
cd ..

mkdir CRR308480_Skin-Het-O-5
cd CRR308480_Skin-Het-O-5
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308480/CRR308480_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308480/CRR308480_r2.fq.gz
cd ..

mkdir CRR308481_Skin-Het-O-6
cd CRR308481_Skin-Het-O-6
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308481/CRR308481_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308481/CRR308481_r2.fq.gz
cd ..

mkdir CRR308482_Skin-Het-O-7
cd CRR308482_Skin-Het-O-7
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308482/CRR308482_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308482/CRR308482_r2.fq.gz
cd ..

mkdir CRR308483_Skin-Het-O-8
cd CRR308483_Skin-Het-O-8
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308483/CRR308483_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR308483/CRR308483_r2.fq.gz
cd ..

<<comment

mkdir Spleen
cd Spleen

mkdir CRR403708_Spleen-Y-CD45.1
cd CRR403708_Spleen-Y-CD45.1
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403708/CRR403708_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403708/CRR403708_r2.fq.gz
cd ..

mkdir CRR403707_Spleen-O-CD45.2
cd CRR403707_Spleen-O-CD45.2
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403707/CRR403707_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403707/CRR403707_r2.fq.gz
cd ..

mkdir CRR403709_Spleen-Y-CD45.2
cd CRR403709_Spleen-Y-CD45.2
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403709/CRR403709_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403709/CRR403709_r2.fq.gz
cd ..

mkdir CRR403706_Spleen-O-CD45.1
cd CRR403706_Spleen-O-CD45.1
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403706/CRR403706_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403706/CRR403706_r2.fq.gz
cd ..

mkdir
cd
wget
wget
cd ..

comment





















 












