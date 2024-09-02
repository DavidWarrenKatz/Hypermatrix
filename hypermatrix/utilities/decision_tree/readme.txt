#create conda enviroment
conda create -n rnaseq_env python=3.8
conda activate rnaseq_env

#install fastdump
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-centos_linux64.tar.gz
tar -vxzf sratoolkit.3.0.10-centos_linux64.tar.gz
export PATH=$PATH:~/sratoolkit.3.0.10-centos_linux64/bin
source ~/.bashrc  # or source ~/.bash_profile
fastq-dump --version

#installl packages needed for RNA pipeline
conda install -c bioconda trimmomatic
conda install -c bioconda tophat
conda install -c bioconda htseq




perl -MNet::FTP -e \
  '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login; $ftp->binary;
  $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
tar -xzvf edirect.tar.gz
rm edirect.tar.gz
builtin echo "export PATH=\$PATH:\$HOME/edirect" >> $HOME/.bash_profile
source $HOME/.bash_profile

esearch -version



#verify installations
trimmomatic -version
tophat --version
htseq-count --version
fastq-dump --version  # For SRA Toolkit






