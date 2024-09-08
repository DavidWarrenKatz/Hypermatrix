[0] set up the conda enrivoment

conda create --name scRNAseq python=3.8
conda activate scRNAseq
conda install -c conda-forge scanpy
pip install scrublet
conda install -c conda-forge jupyterlab
conda install numpy h5py seaborn matplotlib pandas


[1] Both fowards and reverse read FASTQ files for this experiment were downloaded with the following commands. Both files were around 50 GB and this took about 10 minutes. 

wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403693/CRR403693_f1.fq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA004660/CRR403693/CRR403693_r2.fq.gz


[2] Next a matrix was produced by running the following command
There were two files for each experiment:
CRR403693_f1.fq.gz
CRR403693_r2.fq.gz

These files were renamed to work with cell ranger
CRR403693_S1_L001_R2_001.fastq.gz
CRR403693_S1_L001_R1_001.fastq.gz

The following command
cellranger count --id=CRR403693 \
                   --transcriptome=/home/dwk681/workspace/softwareFiles/refdata-gex-GRCh38-and-mm10-2020-A \
                   --fastqs=/home/dwk681/workspace/CRA004660/CRR403693_Liver-Het-O \
                   --sample=CRR403693_Liver-Het-O \
                   --localcores=8 \
                   --localmem=64
by running the runCount.sh file


produced a file CRR403693
which contians these subfiles

-rw-r--r-- 1 dwk681 b1198  741 Oct 21 06:17 _invocation
-rw-r--r-- 1 dwk681 b1198    5 Oct 21 06:17 _jobmode
-rw-r--r-- 1 dwk681 b1198  66K Oct 21 06:17 _mrosource
-rw-r--r-- 1 dwk681 b1198   61 Oct 21 06:17 _versions
-rw-r--r-- 1 dwk681 b1198    2 Oct 21 06:17 _tags
-rw-r--r-- 1 dwk681 b1198   36 Oct 21 06:17 _uuid
drwxrwsr-x 6 dwk681 b1198 4.0K Oct 22 06:29 SC_RNA_COUNTER_CS
-rw-r--r-- 1 dwk681 b1198 1.1M Oct 22 06:30 _vdrkill
drwxrwsr-x 5 dwk681 b1198 4.0K Oct 22 06:30 outs
-rw-r--r-- 1 dwk681 b1198   51 Oct 22 06:30 _timestamp
-rw-r--r-- 1 dwk681 b1198 2.3M Oct 22 06:30 _perf
-rw-r--r-- 1 dwk681 b1198 4.3M Oct 22 06:30 _finalstate
-rw-r--r-- 1 dwk681 b1198 274K Oct 22 06:30 _log
-rw-rw-r-- 1 dwk681 b1198  310 Oct 22 06:30 _cmdline
-rw-rw-r-- 1 dwk681 b1198  17K Oct 22 06:30 _sitecheck
-rw-rw-r-- 1 dwk681 b1198 270K Oct 22 06:30 _filelist
-rw-rw-r-- 1 dwk681 b1198 1.3K Oct 22 06:30 _vdrkill._truncated_
-rw-rw-r-- 1 dwk681 b1198 9.4K Oct 22 06:30 _perf._truncated_
-rw-rw-r-- 1 dwk681 b1198  11M Oct 22 06:31 CRR403693.mri.tgz


