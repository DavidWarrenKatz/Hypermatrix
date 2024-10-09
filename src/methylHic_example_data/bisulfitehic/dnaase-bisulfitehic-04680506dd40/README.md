# Bhmem
A mapping and analysis toolkit for NOMe-HiC, Methyl-HiC, and single-cell Methyl-HiC data


## Pre-requisite:
* JAVA oracle JDK8
* GNU Make 3.81
* gcc 4.8.2
* wget
* maven
* python 3.6
* bismark_genome_preparation in bismark

## Command to install it
```
unzip dnaase-bisulfitehic-5dbb6ce7bc3c.zip
cd dnaase-bisulfitehic-5dbb6ce7bc3c
source ./install.sh
conda create -n bisulfitehic python=3.6
source activate bisulfitehic
pip install numpy
pip install pysam

export JAVA_HOME=/usr/libexec/java_home
```


## Prepare index file for genome and bisulfite converted genomes
* Use bismark_genome_preparation in bismark to prepare the bisulfite converted genomes. It will generate Bisulfite_Genome directory under the same directory as your reference genome, such as mm9.fa.

* Then use bwa index to generate index for each of the converted genomes (Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa and Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa)

* under your reference genome folder, you should have mm9.fa, mm9.fa.fai, mm9.dict, Bisulfite_Genome (directory), Bisulfite_Genome/CT_conversion (directory),  Bisulfite_Genome/GA_conversion (directory), Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa (and its bwa indexed files), and Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa (and its bwa indexed files)


## Analysis 1: Map Methyl-HiC or NOMe-HiC reads (DNA) to the reference genome
```
java -Xmx20G -Djava.library.path=/your/path/to/bisulfitehic/jbwa/src/main/native/ -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem /your/path/to/reference_genome/human_g1k_v37.fa testOut.bam test_R1.fq.gz test_R2.fq.gz -t 5 -rgId test -rgSm sample_id -outputMateDiffChr -buffer 100000
```

### Analysis 1.1: Map single-cell Methyl-HiC reads to the reference genome
```
java -Xmx15G -Djava.library.path=/your/path/to/bisulfitehic/jbwa/src/main/native/ -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/jbwa.jar" main.java.edu.mit.compbio.bisulfitehic.mapping.Bhmem /your/path/to/reference_genome/human_g1k_v37.fa testOut.bam test_R1.fq.gz test_R2.fq.gz -t 1 -rgId test -rgSm sample_id -outputMateDiffChr -buffer 1000 -nonDirectional -pbat -enzymeList hg19_MboI.span_region.bedgraph
```

## Analysis 2: Identify SNPs, DNA methylation, and GpC accessibility in bisulfite sequencing (WGBS, Methyl-HiC, and NOMe-HiC) with high accuracy.
* need to preprocess and index the bam files first
```
samtools sort --threads 5 -T test.name -n $in | samtools fixmate -m --threads 5 - - | samtools sort --threads 5 -T test.cor - | samtools markdup -T test.mdups --threads 5 - - | samtools calmd --threads 5 -b - /your/path/to/reference_genome/human_g1k_v37.fa 2>/dev/null > testOut.calmd.bam
samtools index -@ 10 testOut.calmd.bam
```

* utilize [Bis-SNP.jar](https://github.com/dnaase/Bis-tools/releases/download/Bis-SNP.v0.90/BisSNP-0.90.jar) and [bissnp_easy_usage.pl](https://github.com/dnaase/Bis-tools/blob/master/Bis-SNP/bissnp_easy_usage.pl) from Bis-SNP
* Detail tutorial about Bis-SNP is [here](https://people.csail.mit.edu/dnaase/bissnp2011/)

### Analysis 2.1: Identify DNA methylation, GpC accessibility, and raw SNPs in NOMe-HiC

* Please use 6plus2.strand.bed for the final GCH and HCG methylation analysis
```
perl /your/path/to/Bis-SNP/bissnp_easy_usage.pl --use_bad_mates --nomeseq --mmq 30 --nt 30 --mem 5 /your/path/to/Bis-SNP/BisSNP-0.90.jar testOut.calmd.bam /your/path/to/reference_genome/human_g1k_v37.fa /your/path/to/dbsnp/dbsnp_138.b37.vcf
```

### Analysis 2.2: Identify DNA methylation and raw SNPs in Methyl-HiC
* Please use 6plus2.strand.bed for the final CpG methylation analysis
```
perl /your/path/to/Bis-SNP/bissnp_easy_usage.pl --use_bad_mates --mmq 30 --nt 30 --mem 5 /your/path/to/Bis-SNP/BisSNP-0.90.jar testOut.calmd.bam /your/path/to/reference_genome/human_g1k_v37.fa /your/path/to/dbsnp/dbsnp_138.b37.vcf
```

### Analysis 2.3: Identify SNPs with comparable high accuracy as that in WGS
* Please use snp.filtered.sort.vcf from Analysis 2.1 or Analysis 2.2 for the further SNP calls
1. bgzip and index vcf file
```
bgzip -c testOut.calmd.snp.filtered.sort.vcf > testOut.calmd.snp.filtered.sort.vcf.gz
tabix -p vcf -b 2 -e 2 testOut.calmd.snp.filtered.sort.vcf.gz
```
2. only use SNP loci for the analysis and split the vcf file into single chromosome
```
perl -e '@cs=(1..22);foreach $c(@cs){`bcftools view -v snps -g ^miss -O z -o testOut.bissnp_filter.snp_only.chr$c.vcf.gz -r $c testOut.calmd.snp.filtered.sort.vcf.gz && tabix -p vcf -b 2 -e 2 testOut.bissnp_filter.snp_only.chr$c.vcf.gz`;print "$c\n";}'
```
3. prepare SNP chunk. need to download [sites.vcf](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) files from 1000 genome project.
```
perl -e '%chr_size=();open(F,"</your/path/to/reference_genome/human_g1k_v37.autosome.chrom.sizes");while(<F>){chomp;@f=split "\t";$chr_size{$f[0]}=$f[1];}close(F);@cs=(1..22);foreach my $c(@cs){print "$c\n";open(O, ">chunk.testOut.bissnp_filter.1000GP.EUR.noNA12878.chr$c.txt") or die;for($s=0;$s<$chr_size{$c};$s+=20000000){$e=$chr_size{$c} >= ($s+20000000) ? ($s+20000000) : $chr_size{$c}; $num=`bcftools view /your/path/to/1kg_panel/1000GP.EUR.chr$c.noNA12878.sites.vcf.gz -r ${c}:${s}-${e} | bedtools intersect -wa -u -a testOut.bissnp_filter.snp_only.chr$c.vcf.gz -b - | wc -l`;chomp($num);$total=`bcftools view -H testOut.bissnp_filter.snp_only.chr$c.vcf.gz -r ${c}:${s}-${e} | wc -l`;chomp($total);if($num>3 and $num/$total>0.5){print O "$c\t$s\t$e\t$num\t$total\n";}}close(O);}'
```
4. perform imputation and phasing of SNPs at each of the SNP chunk. need to download [genotypes.vcf](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) files from 1000 genome project.
* nedd to download and install [eagle](https://alkesgroup.broadinstitute.org/Eagle/), [minimac4](https://github.com/statgen/Minimac4) and their related files
```
perl -e '@cs=(1..22);foreach my $c(@cs){open(F,"<chunk.testOut.bissnp_filter.1000GP.EUR.noNA12878.chr$c.txt");$i=0;while(<F>){chomp;my @f=split "\t";$s=$f[1]+1;$e=$f[2];print "$i\t$s\t$e\n";`/your/path/to/Eagle/default/eagle --numThreads 10 --chrom $c --bpStart $s --bpEnd $e --bpFlanking 5000000 --allowRefAltSwap --geneticMapFile /your/path/to/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --outPrefix testOut.bissnp_filter.eagle2_prephased.chr$c.$i --vcfRef /your/path/to/1kg_panel/1000GP.EUR.chr$c.noNA12878.bcf --vcfTarget testOut.bissnp_filter.snp_only.chr$c.vcf.gz --vcfOutFormat z && /your/path/to/Minimac4/bin/minimac4 --refHaps /your/path/to/1kg_panel/refPanel.1000GP.EUR.chr$c.noNA12878.m3vcf.gz --haps testOut.bissnp_filter.eagle2_prephased.chr$c.$i.vcf.gz --mapFile /your/path/to/Eagle/default/tables/genetic_map_hg19_withX.txt.gz --prefix testOut.bissnp_filter.eagle2_prephased.minimac4_impute.chr$c.$i --cpus 10 --window 500000 --noPhoneHome --format GT,DS,GP --allTypedSites --meta --minRatio 0.00001 --chr $c --start $s --end $e`;$i++;}close(F);}'
```
5. merge all SNP files after imputation and phasing.
```
bcftools concat -O z testOut.bissnp_filter.eagle2_prephased.minimac4_impute.chr*.*.dose.vcf.gz | bcftools sort -m 7680M -O b -o testOut.bissnp_filter.eagle2_prephased.minimac4_impute.autosomes.bcf - 
bcftools index testOut.bissnp_filter.eagle2_prephased.minimac4_impute.autosomes.bcf
bcftools view -m 2 -M 2 -v snps -i 'R2>0.3 & MAF>0.01' -O b -o testOut.bissnp_filter.eagle2_prephased.minimac4_impute.filtered.autosomes.bcf testOut.bissnp_filter.eagle2_prephased.minimac4_impute.autosomes.bcf
bcftools index testOut.bissnp_filter.eagle2_prephased.minimac4_impute.filtered.autosomes.bcf
```
6. check the concordance with vcf file from 1000G's WGS result by GATK4.
```
gatk Concordance -R /your/path/to/reference_genome/human_g1k_v37.fa -eval testOut.bissnp_filter.eagle2_prephased.minimac4_impute.filtered.autosomes.vcf.gz --truth NA12878.eagle2_prephased.minimac4_impute.filtered.autosomes.vcf.gz --summary testOut_vs_1000G_WGS.bissnp_filter.impute.summary.tsv 
```

## Analysis 3: Characterize 3D genome information from the mapped bam files
* convert bam files to .hic file
```
python /your/path/to/bisulfitehic/src/python/sam2juicer.py -s testOut.calmd.bam -f /your/path/to/juicer/restriction_sites/hg19_DpnII.txt | gzip -fc > testOut.hic.txt.gz
java -Xmx20G -jar /your/path/to/juicer/scripts/common/juicer_tools.jar pre -d -q 30 -f /your/path/to/juicer/restriction_sites/hg19_DpnII.txt testOut.hic.txt.gz testOut.mapQ30.hic hg19
```

* call TADs from .hic file by [TopDom](https://cran.r-project.org/web/packages/TopDom/index.html)
```
Rscript /your/path/to/bisulfitehic/src/R/call_topdom_on_hic.R testOut.mapQ30.hic testOut.mapQ30.topdom_default.tsv
```
* call chromatin loops from .hic file by [Mustache](https://github.com/ay-lab/mustache)
```
conda activate mustache
python /your/path/to/mustache/mustache.py -f testOut.mapQ30.hic -ch 22 -p 5 -r 25000 -cz /your/path/to/reference_genome/human_g1k_v37.chrom.sizes -norm KR -d 5000000 -o ./testOut.chr22.mustache.tsv
```
* call compartment from .hic file by [FAN-C](https://vaquerizaslab.github.io/fanc/fanc-executable/fanc-analyse-hic/ab_compartments.html) or [Juicer](https://github.com/aidenlab/juicer/wiki/Eigenvector)

## Analysis 4: Identify the long-range epigenetic concordance at specific regions

### Analysis 4.1: long-range (>20kb) GCH methylation correlation in NOMe-HiC
```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups GCH_concordance.testOut.mustache testOut.all_chr.mustache.bedpe testOut.b37.calmd.bam -methyPatternSearch GCH -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate -useGeneralBedPe
```
### Analysis 4.2: long-range (>20kb) WCG methylation correlation in NOMe-HiC
```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups WCG_concordance.testOut.mustache testOut.all_chr.mustache.bedpe testOut.b37.calmd.bam -methyPatternSearch WCG -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate -useGeneralBedPe
```

### Analysis 4.3: short-range (<1lb) GCH methylation correlation in NOMe-HiC
```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups GCH_concordance_short.testOut.mustache testOut.all_chr.mustache.bedpe testOut.b37.calmd.bam -methyPatternSearch GCH -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate -onlyShort 1000 -randomSampleReads 100 -useGeneralBedPe
```
### Analysis 4.4: long-range (>20kb) GCH methylation correlation but in shuffled reads in NOMe-HiC
```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.MethyCorAcrossHiccups GCH_concordance_shuffle.testOut.mustache testOut.all_chr.mustache.bedpe testOut.b37.calmd.bam -methyPatternSearch GCH -methyPatternCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate -randomHiccupsInput -randomSampleReads 100 -useGeneralBedPe
```

* GCH_concordance.testOut.mustache.methyCor.txt.gz can used as an input to calculate the correlation level in R
* Example R script (need to install [psych](https://cran.r-project.org/web/packages/psych/index.html) package):
```
library(psych)

minCT=3
d<-read.table("GCH_concordance_short.testOut.mustache.methyCor.txt.gz",sep="\t",header=F)
cor.sameReads.short<-cor.test(d[d[,13]>=minCT&d[,15]>=minCT,12],d[d[,13]>=minCT&d[,15]>=minCT,14])$estimate
s<-cbind(d[d[,13]>=minCT&d[,15]>=minCT,12],d[d[,13]>=minCT&d[,15]>=minCT,14])
s[s<0.5]=0
s[s>=0.5]=1
data<-matrix(c(length(s[s[,1]==0&s[,2]==0,2]), length(s[s[,1]==0&s[,2]==1,2]), length(s[s[,1]==1&s[,2]==0,2]), length(s[s[,1]==1&s[,2]==1,2])), nrow = 2) 
phi.sameReads.short<-phi(data, digits = 3)
t_cor.sameReads.short<-tetrachoric(data)[1]$rho


d<-read.table("GCH_concordance.testOut.mustache.methyCor.txt.gz",sep="\t",header=F)
cor.sameReads.long<-cor.test(d[d[,13]>=minCT&d[,15]>=minCT,12],d[d[,13]>=minCT&d[,15]>=minCT,14])$estimate
s<-cbind(d[d[,13]>=minCT&d[,15]>=minCT,12],d[d[,13]>=minCT&d[,15]>=minCT,14])
s[s<0.5]=0
s[s>=0.5]=1
data<-matrix(c(length(s[s[,1]==0&s[,2]==0,2]), length(s[s[,1]==0&s[,2]==1,2]), length(s[s[,1]==1&s[,2]==0,2]), length(s[s[,1]==1&s[,2]==1,2])), nrow = 2) 
phi.sameReads.long<-phi(data, digits = 3)
t_cor.sameReads.long<-tetrachoric(data)[1]$rho

d<-read.table("GCH_concordance_shuffled.testOut.mustache.methyCor.txt.gz",sep="\t",header=F)
cor.acrossReads<-cor.test(d[d[,13]>=minCT&d[,15]>=minCT,12],d[d[,13]>=minCT&d[,15]>=minCT,14])$estimate
s<-cbind(d[d[,13]>=minCT&d[,15]>=minCT,12],d[d[,13]>=minCT&d[,15]>=minCT,14])
s[s<0.5]=0
s[s>=0.5]=1
data<-matrix(c(length(s[s[,1]==0&s[,2]==0,2]), length(s[s[,1]==0&s[,2]==1,2]), length(s[s[,1]==1&s[,2]==0,2]), length(s[s[,1]==1&s[,2]==1,2])), nrow = 2) 
phi.acrossReads<-phi(data, digits = 3)
t_cor.acrossReads<-tetrachoric(data)[1]$rho


data<-cbind(c(cor.sameReads.short, cor.sameReads.long, cor.acrossReads),c(phi.sameReads.short, phi.sameReads.long, phi.acrossReads),c(t_cor.sameReads.short, t_cor.sameReads.long, t_cor.acrossReads))
rownames(data)<-c("Short (<1kb)","Long (>20kb)","Shuffle (>20kb)")
colnames(data)<-c("Pearson","Phi","Tetrachoric")
pdf("NOMeHiC.testOut.mustache_25kb_res.GCH_pearson_phi_tet_cor.pdf", paper="special", height=5, width=5)
par(mar = c(4,4,3,1))
barplot(height = data, beside = TRUE, ylab="Correlation Coefficient",ylim=c(0,1),xlab="",main="",col=c("black","orange","blue"))
legend("topleft",c("Local (<1kb)","Long-range (>20kb)","Shuffle (>20kb)"),fill=c("black","orange","blue"),cex=0.75)
dev.off()
```

## Analysis 5: Identify the long-range allele-specific epigenetic events (NOMe-HiC) at chr22
```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.LongRangeAsm /your/path/to/phased_vcf_files/testOut.vcf.gz testOut.b37.calmd.bam testOut.ASM_GCH_HCG_20kb_minCov2.chr22.tab.gz testOut.ASM_GCH_HCG_20kb_minCov2.ref.chr22.tab.gz testOut.ASM_GCH_HCG_20kb_minCov2.alt.chr22.tab.gz -minDist 20000 -sample 0 -coverageRef 2 -coverageAlt 2 -methyPatternSearch GCH -methyPatternCytPos 2 -methyPattern2ndSearch HCG -methyPattern2ndCytPos 2 -bsConvPatternSearch WCH -bsConvPatternCytPos 2 -useBadMate  -minMapQ 30 -region 22
```
* testOut.ASM_GCH_HCG_20kb_minCov2.chr22.tab.gz is the output file, 
* Merged files from all the chromosome:
```
zcat testOut.ASM_GCH_HCG_20kb_minCov2.chr*.tab.gz | grep -v "#" | gzip -c > testOut.ASM_GCH_HCG_20kb_minCov2.all_chr.tab.gz
```
* testOut.ASM_GCH_HCG_20kb_minCov2.all_chr.tab.gz can used as an input to further identify the significant allele-specific GCH methylated regions in R
* Example R script:
```
x<-read.table("testOut.ASM_GCH_HCG_20kb_minCov2.all_chr.tab.gz",sep="\t",header=F)
rownames(x)<-paste(x[,1],x[,2],sep=":")
x.meth<-as.matrix(x[,c(8,13,17,20,25,29)]) #only consider GCH first
x<-x[rowSums(is.na(x.meth))==0 & rowSums(is.infinite(x.meth))==0 & rowSums(x.meth<0)==0,]

#AlleleA, AlleleB, AlleleA_end2, AlleleB_end2, AlleleA_local, AlleleB_local, AlleleA_local_HCG, AlleleB_local_HCG, AlleleA_HCG, AlleleB_HCG, AlleleA_end2_HCG, Alleleb_end2_HCG
d<-cbind(as.numeric(x[,6]),as.numeric(x[,7]),#AlleleA
as.numeric(x[,11]),as.numeric(x[,12]),#AlleleB
as.numeric(x[,15]),as.numeric(x[,16]),#AlleleA_end2
as.numeric(x[,18]),as.numeric(x[,19]),#AlleleB_end2
as.numeric(x[,23]),as.numeric(x[,24]),#AlleleA_local
as.numeric(x[,27]),as.numeric(x[,28]),#AlleleB_local
as.numeric(x[,30]),as.numeric(x[,31]),#AlleleA_local_HCG
as.numeric(x[,33]),as.numeric(x[,34]),#AlleleB_local_HCG
as.numeric(x[,36]),as.numeric(x[,37]),#AlleleA_HCG
as.numeric(x[,39]),as.numeric(x[,40]),#AlleleB_HCG
as.numeric(x[,42]),as.numeric(x[,43]),#AlleleA_end2_HCG
as.numeric(x[,45]),as.numeric(x[,46]) #Alleleb_end2_HCG
)
rownames(d)<-paste(x[,1],x[,2],sep=":")

getScore<-function(s){
c(fisher.test(matrix(c(s[1],s[2]-s[1],s[3],s[4]-s[3]),nrow=2))$p.value, #end1
fisher.test(matrix(c(s[5],s[6]-s[5],s[7],s[8]-s[7]),nrow=2))$p.value,  #end2
fisher.test(matrix(c(s[5]+s[1],s[6]-s[5]+s[2]-s[1],s[7]+s[3],s[8]-s[7]+s[4]-s[3]),nrow=2))$p.value,   #end1+end2
fisher.test(matrix(c(s[9],s[10]-s[9],s[11],s[12]-s[11]),nrow=2))$p.value, #local
abs(s[1]/s[2]-s[3]/s[4]),
abs(s[5]/s[6]-s[7]/s[8]),
abs((s[5]+s[1])/(s[6]+s[2])-(s[7]+s[3])/(s[4]+s[8]))
)
}
test.score<-t(apply(d,1,getScore))

test.qvalue<-cbind(p.adjust(test.score[,1],method="BH"),
p.adjust(test.score[,2],method="BH"),
p.adjust(test.score[,3],method="BH"),
p.adjust(test.score[,4],method="BH")
)

#FDR<0.05
ASM.1stEndAnd2ndEnd<-test.qvalue[test.qvalue[,1]<0.05 & test.qvalue[,2]<0.05 & test.qvalue[,4]<0.05,]
ASM.1stEndOnly<-test.qvalue[test.qvalue[,1]<0.05 & test.qvalue[,2]>0.95 & test.qvalue[,4]<0.05,]
ASM.2ndEndOnly<-test.qvalue[test.qvalue[,2]<0.05 & test.qvalue[,1]>0.95 & test.qvalue[,4]>0.95,]
ASM.1stEnd<-test.qvalue[test.qvalue[,4]<0.05,]

write.table(x[rownames(ASM.1stEndAnd2ndEnd),],"GM12878_bT_WGBS_merged_final.ASM_GCH_HCG_20kb_minCov2.merge_all.ASM.1stEndAnd2ndEnd.and_local.fdr005.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(x[rownames(ASM.1stEndOnly),],"GM12878_bT_WGBS_merged_final.ASM_GCH_HCG_20kb_minCov2.merge_all.ASM.1stEndOnly.and_local.fdr005.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(x[rownames(ASM.2ndEndOnly),],"GM12878_bT_WGBS_merged_final.ASM_GCH_HCG_20kb_minCov2.merge_all.ASM.2ndEndOnly.and_local.fdr005.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(x[rownames(ASM.1stEnd),],"GM12878_bT_WGBS_merged_final.ASM_GCH_HCG_20kb_minCov2.merge_all.ASM.1stEnd.and_local.fdr005.txt",sep="\t",quote=F,row.names=F,col.names=F)
```

## Help information about specific options in each command can be obtained by simply type: 
```
java -Xmx20G -cp "/your/path/to/bisulfitehic/target/bisulfitehic-0.38-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/ccInference-0.19-jar-with-dependencies.jar:/your/path/to/bisulfitehic/lib/gatk-package-distribution-3.3.jar" main.java.edu.mit.compbio.bisulfitehic.utils.LongRangeAsm
```

## Citations:
### For NOMe-HiC application
Fu H, Zheng H, Muglia LJ, Wang L, Liu Y. GTAGMe-seq: joint profiling of genetic variants, DNA methylation, GpC methyltransferase footprints, and 3D genome in the same DNA molecules. BioRxiv preprint. 2022 Mar; doi: https://doi.org/10.1101/2022.03.29.486102

### For Methyl-HiC application
Li G*, Liu Y*, Zhang Y, Kubo N, Yu M, Fang R, Kellis M, Ren B. Joint profiling of DNA methylation and chromatin architecture in single cells. Nat Methods. 2019 Oct;16(10):991-993. doi: https://doi.org/10.1038/s41592-019-0502-z. PubMed PMID: 31384045; PubMed Central PMCID: PMC6765429






