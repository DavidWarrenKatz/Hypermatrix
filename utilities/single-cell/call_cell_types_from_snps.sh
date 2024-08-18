#!/bin/bash
# Author: "Author Name"
# File: call_cell_types_from_snps.sh
# Date: 08/18/2024
# Description: Script for genomic data analysis: merging, SNP calling, and imputation


# Variables and filepaths
author="Author Name"
file_description="Genomic Data Analysis"
date=$(date +%Y-%m-%d)
project_directory="/jet/home/dnaase/startup/projects/hexa_seq/scNOMeHiC_RNA_IMR90_GM_20210701"
bam_directory="$project_directory/bam/DNA/merged_cluster"
genome_reference="/jet/home/dnaase/shared/data/genomes/hg19/human_g1k_v37.fa"
dbsnp_reference="/jet/home/dnaase/shared/data/dbsnp/dbsnp_138.b37.vcf"
software_directory="/jet/home/dnaase/startup/software"
bis_tools_directory="$software_directory/bistools/Bis-tools/Bis-SNP"
genetic_map="$software_directory/Eagle/default/tables/genetic_map_hg19_withX.txt.gz"
chunk_size=20000000  # chunk size for imputation

# Functions 
# --------------------------------------------------------- 
# Function to create bam list
function create_bam_list {
    color=$1
    grep "$color" name.order.HCG_methy.with_color.txt | cut -f2 | perl -ne 'chomp; $in=$_; $in.=".b37.calmd.bam"; print "$in\n";' > "${color}_cluster.bam_List.txt"
}

# Function to merge_and_index_bam using bamtools and samtools
# Consensus bam file 
function merge_and_index_bam {
    color=$1
    list_file="${color}_cluster.bam_List.txt"
    output_bam="merge_${color}_cluster.b37.calmd.bam"
    perl ~/ssubexc.pl "bamtools merge -list $list_file -out $output_bam && samtools index -@ 5 $output_bam" "${color}_cluster" 10 5
}

# Function to call SNPs using Bis-SNP, designed for bisulfite sequencing data
# Variants output vcf 

function call_snps {
    bam_file=$1
    output_prefix=$(basename "$bam_file" .bam)
    perl ~/ssubexc.pl "perl $bis_tools_directory/bissnp_easy_usage.test.pl --use_bad_mates --nomeseq --mmq 30 --nt 6 --mem 40 $bis_tools_directory/BisSNP-0.90.jar $bam_file $genome_reference $dbsnp_reference" "bissnp.map30.$output_prefix" 40 20
}



function split_and_impute_vcf {
    color=$1
    vcf="merge_${color}_cluster.b37.calmd.snp.filtered.sort.vcf"
    bgzip -c "$vcf" > "${vcf}.gz"
    tabix -p vcf -b2 -e2 "${vcf}.gz"
    for c in {1..22}
    do
        output_vcf="merge_${color}_cluster.b37.bis_snp.filtered.chr$c.vcf.gz"
        bcftools view -v snps -O z -o "$output_vcf" -r $c "${vcf}.gz"
        bcftools index "$output_vcf"
    done
}


# MAIN 
# --------------------------------------------------------------
create_bam_list "red"
create_bam_list "blue"
merge_and_index_bam "red"
merge_and_index_bam "blue"

ls merge_*.calmd.bam | while read bam_file; do 
    call_snps "$bam_files"
done 

# gatk 
split_and_impute_vcf "red"
split_and_impute_vcf "blue"


# # OLD SCRIPT VERSION TO BE MOVED TO ARCHIVE 
# #merge_red: 74340625 high quality reads (~3X), merge_blue: 95972988 high quality (~4X).
# ##use merged bam file in each cluster, and then call SNPs from each cluster and identify their concordance with each genome
# #/jet/home/dnaase/startup/projects/hexa_seq/scNOMeHiC_RNA_IMR90_GM_20210701/bam/DNA/merged_cluster
# #  1. takes methylation 
# # clusters methylation 
# # 2 groups 
# grep red name.order.HCG_methy.with_color.txt | cut -f2 | perl -ne 'chomp;$in=$_;$in.=".b37.calmd.bam";print "$in\n";' > red_cluster.bam_List.txt
# grep blue name.order.HCG_methy.with_color.txt | cut -f2 | perl -ne 'chomp;$in=$_;$in.=".b37.calmd.bam";print "$in\n";' > blue_cluster.bam_List.txt

# perl ~/ssubexc.pl "bamtools merge -list red_cluster.bam_List.txt -out merge_red_cluster.b37.calmd.bam && samtools index -@ 5 merge_red_cluster.b37.calmd.bam" merge_red_cluster 10 5
# perl ~/ssubexc.pl "bamtools merge -list blue_cluster.bam_List.txt -out merge_blue_cluster.b37.calmd.bam && samtools index -@ 5 merge_blue_cluster.b37.calmd.bam" merge_blue_cluster 10 5

# ##call SNPs from the file
# STEP2 : SNP CALLING 
# ls merge_*.calmd.bam | perl -ne 'chomp;$in=$_;`perl ~/ssubexc.pl "perl /jet/home/dnaase/startup/software/bistools/Bis-tools/Bis-SNP/bissnp_easy_usage.test.pl --use_bad_mates --nomeseq --mmq 30 --nt 6 --mem 40 /jet/home/dnaase/startup/software/bistools/Bis-tools/Bis-SNP/BisSNP-0.90.jar $in /jet/home/dnaase/shared/data/genomes/hg19/human_g1k_v37.fa /jet/home/dnaase/shared/data/dbsnp/dbsnp_138.b37.vcf" bissnp.map30.$in 40 20`;'


# ##split vcf to 22 chrs and impute it by Minimac4 with noNA12878 vcf panel...
##

# perl -e '@cs=(1..22);foreach $c(@cs){`bgzip -c merge_blue_cluster.b37.calmd.snp.filtered.sort.vcf > merge_blue_cluster.b37.calmd.snp.filtered.sort.vcf.gz && tabix -p vcf -b2 -e2 merge_blue_cluster.b37.calmd.snp.filtered.sort.vcf.gz && bcftools view -v snps -O z -o merge_blue_cluster.b37.bis_snp.filtered.chr$c.vcf.gz -r $c merge_blue_cluster.b37.calmd.snp.filtered.sort.vcf.gz && bcftools index merge_blue_cluster.b37.bis_snp.filtered.chr$c.vcf.gz`;print "$c\n";}'
# perl -e '@cs=(1..22);foreach $c(@cs){`bgzip -c merge_red_cluster.b37.calmd.snp.filtered.sort.vcf > merge_red_cluster.b37.calmd.snp.filtered.sort.vcf.gz && tabix -p vcf -b2 -e2 merge_red_cluster.b37.calmd.snp.filtered.sort.vcf.gz && bcftools view -v snps -O z -o merge_red_cluster.b37.bis_snp.filtered.chr$c.vcf.gz -r $c merge_red_cluster.b37.calmd.snp.filtered.sort.vcf.gz && bcftools index merge_red_cluster.b37.bis_snp.filtered.chr$c.vcf.gz`;print "$c\n";}'


# ##prepare the chuncks in each chrs (at least 3 variants and >20% variants overlapped with ref in each chunk).
# # STEP 4 BREAK CHROMSOME INTO MANAGEABLE CHUNKS CALCULATE OVERLAPS WITH REFERENCE , TAKE OUT PUT STATICS FOR POSSIBLE IMPUTATION 
# perl -e '%chr_size=();open(F,"</jet/home/dnaase/shared/data/genomes/hg19/human_g1k_v37.autosome.chrom.sizes");while(<F>){chomp;@f=split "\t";$chr_size{$f[0]}=$f[1];}close(F);@cs=(1..22);foreach my $c(@cs){open(O, ">chunk.merge_red_cluster.1000GP.EUR.noNA12878.chr$c.txt") or die;for($s=0;$s<$chr_size{$c};$s+=20000000){$e=$chr_size{$c} >= ($s+20000000) ? ($s+20000000) : $chr_size{$c}; $num=`bcftools view /jet/home/dnaase/shared/data/genotype/1kg_panel/1000GP.EUR.chr$c.noNA12878.sites.vcf.gz -r ${c}:${s}-${e} | bedtools intersect -wa -u -a merge_red_cluster.b37.bis_snp.filtered.chr$c.vcf.gz -b - | wc -l`;chomp($num);$total=`bcftools view -H merge_red_cluster.b37.bis_snp.filtered.chr$c.vcf.gz -r ${c}:${s}-${e} | wc -l`;chomp($total);if($num>3 and $num/$total>0.2){print O "$c\t$s\t$e\t$num\t$total\n";}}close(O);}'
# perl -e '%chr_size=();open(F,"</jet/home/dnaase/shared/data/genomes/hg19/human_g1k_v37.autosome.chrom.sizes");while(<F>){chomp;@f=split "\t";$chr_size{$f[0]}=$f[1];}close(F);@cs=(1..22);foreach my $c(@cs){open(O, ">chunk.merge_blue_cluster.1000GP.EUR.noNA12878.chr$c.txt") or die;for($s=0;$s<$chr_size{$c};$s+=20000000){$e=$chr_size{$c} >= ($s+20000000) ? ($s+20000000) : $chr_size{$c}; $num=`bcftools view /jet/home/dnaase/shared/data/genotype/1kg_panel/1000GP.EUR.chr$c.noNA12878.sites.vcf.gz -r ${c}:${s}-${e} | bedtools intersect -wa -u -a merge_blue_cluster.b37.bis_snp.filtered.chr$c.vcf.gz -b - | wc -l`;chomp($num);$total=`bcftools view -H merge_blue_cluster.b37.bis_snp.filtered.chr$c.vcf.gz -r ${c}:${s}-${e} | wc -l`;chomp($total);if($num>3 and $num/$total>0.2){print O "$c\t$s\t$e\t$num\t$total\n";}}close(O);}'

# #perl -e '@cs=(1..22);foreach my $c(@cs){open(F,"<chunk.merge_red_cluster.1000GP.EUR.noNA12878.chr$c.txt");$i=0;while(<F>){chomp;my @f=split "\t";$s=$f[1]+1;$e=$f[2];print "$i\t$s\t$e\n";`/jet/home/dnaase/startup/software/Eagle/default/eagle --numThreads 10 --chrom $c --bpStart $s --bpEnd $e --bpFlanking 5000000 --allowRefAltSwap --geneticMapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --outPrefix merge_red_cluster.eagle2_prephased.chr$c.$i --vcfRef /jet/home/dnaase/shared/data/genotype/1kg_panel/1000GP.EUR.chr$c.noNA12878.bcf --vcfTarget merge_red_cluster.b37.bis_snp.filtered.chr$c.vcf.gz --vcfOutFormat z && /jet/home/dnaase/startup/software/Minimac4/bin/minimac4 --refHaps /jet/home/dnaase/shared/data/genotype/1kg_panel/refPanel.1000GP.EUR.chr$c.noNA12878.m3vcf.gz --haps merge_red_cluster.eagle2_prephased.chr$c.$i.vcf.gz --mapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --prefix merge_red_cluster.eagle2_prephased.minimac4_impute.chr$c.$i --cpus 10 --window 500000 --noPhoneHome --format GT,DS,GP --allTypedSites --meta --minRatio 0.00001 --chr $c --start $s --end $e`;$i++;}close(F);}'
# #perl -e '@cs=(1..22);foreach my $c(@cs){open(F,"<chunk.merge_blue_cluster.1000GP.EUR.noNA12878.chr$c.txt");$i=0;while(<F>){chomp;my @f=split "\t";$s=$f[1]+1;$e=$f[2];print "$i\t$s\t$e\n";`/jet/home/dnaase/startup/software/Eagle/default/eagle --numThreads 10 --chrom $c --bpStart $s --bpEnd $e --bpFlanking 5000000 --allowRefAltSwap --geneticMapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --outPrefix merge_blue_cluster.eagle2_prephased.chr$c.$i --vcfRef /jet/home/dnaase/shared/data/genotype/1kg_panel/1000GP.EUR.chr$c.noNA12878.bcf --vcfTarget merge_blue_cluster.b37.bis_snp.filtered.chr$c.vcf.gz --vcfOutFormat z && /jet/home/dnaase/startup/software/Minimac4/bin/minimac4 --refHaps /jet/home/dnaase/shared/data/genotype/1kg_panel/refPanel.1000GP.EUR.chr$c.noNA12878.m3vcf.gz --haps merge_blue_cluster.eagle2_prephased.chr$c.$i.vcf.gz --mapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --prefix merge_blue_cluster.eagle2_prephased.minimac4_impute.chr$c.$i --cpus 10 --window 500000 --noPhoneHome --format GT,DS,GP --allTypedSites --meta --minRatio 0.00001 --chr $c --start $s --end $e`;$i++;}close(F);}'

# perl -e '@cs=(1..22);foreach my $c(@cs){open(F,"<chunk.merge_red_cluster.1000GP.EUR.noNA12878.chr$c.txt");$i=0;while(<F>){chomp;my @f=split "\t";$s=$f[1]+1;$e=$f[2];print "$i\t$s\t$e\n";`perl /jet/home/dnaase/ssubexc.pl "/jet/home/dnaase/startup/software/Eagle/default/eagle --numThreads 10 --chrom $c --bpStart $s --bpEnd $e --bpFlanking 5000000 --allowRefAltSwap --geneticMapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --outPrefix merge_red_cluster.eagle2_prephased.chr$c.$i --vcfRef /jet/home/dnaase/shared/data/genotype/1kg_panel/1000GP.EUR.chr$c.noNA12878.bcf --vcfTarget merge_red_cluster.b37.bis_snp.filtered.chr$c.vcf.gz --vcfOutFormat z && /jet/home/dnaase/startup/software/Minimac4/bin/minimac4 --refHaps /jet/home/dnaase/shared/data/genotype/1kg_panel/refPanel.1000GP.EUR.chr$c.noNA12878.m3vcf.gz --haps merge_red_cluster.eagle2_prephased.chr$c.$i.vcf.gz --mapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --prefix merge_red_cluster.eagle2_prephased.minimac4_impute.chr$c.$i --cpus 10 --window 500000 --noPhoneHome --format GT,DS,GP --allTypedSites --meta --minRatio 0.00001 --chr $c --start $s --end $e" merge_red_impute_${c}_$i 20 10`;$i++;}close(F);}'
# perl -e '@cs=(1..22);foreach my $c(@cs){open(F,"<chunk.merge_blue_cluster.1000GP.EUR.noNA12878.chr$c.txt");$i=0;while(<F>){chomp;my @f=split "\t";$s=$f[1]+1;$e=$f[2];print "$i\t$s\t$e\n";`perl /jet/home/dnaase/ssubexc.pl "/jet/home/dnaase/startup/software/Eagle/default/eagle --numThreads 10 --chrom $c --bpStart $s --bpEnd $e --bpFlanking 5000000 --allowRefAltSwap --geneticMapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --outPrefix merge_blue_cluster.eagle2_prephased.chr$c.$i --vcfRef /jet/home/dnaase/shared/data/genotype/1kg_panel/1000GP.EUR.chr$c.noNA12878.bcf --vcfTarget merge_blue_cluster.b37.bis_snp.filtered.chr$c.vcf.gz --vcfOutFormat z && /jet/home/dnaase/startup/software/Minimac4/bin/minimac4 --refHaps /jet/home/dnaase/shared/data/genotype/1kg_panel/refPanel.1000GP.EUR.chr$c.noNA12878.m3vcf.gz --haps merge_blue_cluster.eagle2_prephased.chr$c.$i.vcf.gz --mapFile /jet/home/dnaase/startup/software/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --prefix merge_blue_cluster.eagle2_prephased.minimac4_impute.chr$c.$i --cpus 10 --window 500000 --noPhoneHome --format GT,DS,GP --allTypedSites --meta --minRatio 0.00001 --chr $c --start $s --end $e" merge_blue_impute_${c}_$i 20 10`;$i++;}close(F);}'


# ##merge imputed vcf file
# bcftools concat -O u merge_red_cluster.eagle2_prephased.minimac4_impute.chr*.*.dose.vcf.gz | bcftools sort -m 7680M -O u - -| bcftools view -m 2 -M 2 -v snps -i 'R2>0.3 & MAF>0.01' -O z -o merge_red_cluster.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz - && tabix -p vcf -b2 -e2 merge_red_cluster.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz
# bcftools concat -O u merge_blue_cluster.eagle2_prephased.minimac4_impute.chr*.*.dose.vcf.gz | bcftools sort -m 7680M -O u - -| bcftools view -m 2 -M 2 -v snps -i 'R2>0.3 & MAF>0.01' -O z -o merge_blue_cluster.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz - && tabix -p vcf -b2 -e2 merge_blue_cluster.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz

# #bcftools view -m 2 -M 2 -v snps -i 'R2>0.7 & MAF>0.01' -O z -o merge_red_cluster.eagle2_prephased.minimac4_impute.filtered.chr22.vcf.gz merge_red_cluster.eagle2_prephased.minimac4_impute.chr22.vcf.gz && tabix -p vcf -b2 -e2 merge_red_cluster.eagle2_prephased.minimac4_impute.filtered.chr22.vcf.gz
# #bcftools view -m 2 -M 2 -v snps -i 'R2>0.7 & MAF>0.01' -O z -o merge_blue_cluster.eagle2_prephased.minimac4_impute.filtered.chr22.vcf.gz merge_blue_cluster.eagle2_prephased.minimac4_impute.chr22.vcf.gz && tabix -p vcf -b2 -e2 merge_blue_cluster.eagle2_prephased.minimac4_impute.filtered.chr22.vcf.gz



# #should use gatk to evaluate
# gatk Concordance -R /jet/home/dnaase/shared/data/genomes/hg19/human_g1k_v37.fa -eval merge_red_cluster.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz --truth /jet/home/dnaase/shared/data/genotype/IMR90/illumina/IMR90.WGS_2D.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz --summary merge_red_cluster.IMR90_WGS_2D.summary.tsv 
# gatk Concordance -R /jet/home/dnaase/shared/data/genomes/hg19/human_g1k_v37.fa -eval merge_blue_cluster.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz --truth /jet/home/dnaase/shared/data/genotype/IMR90/illumina/IMR90.WGS_2D.eagle2_prephased.minimac4_impute.filtered.autosome.vcf.gz --summary merge_blue_cluster.IMR90_WGS_2D.summary.tsv 


