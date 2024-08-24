        conda create --name multiomics10
        conda activate multiomics10
        pip install hic-straw
        conda install nomkl
        conda install numpy
        conda install matplotlib
        conda install seaborn
        conda install requests
        conda install h5py
        conda install scipy
        pip install pyBigWig
        pip install fanc




#STEP1:	CREATE AND SETUP CONDA ENVIRONMENT
#I first created the conda enviroment multiomics6 with the following commands

	module load python/anaconda3.6
	conda create -n multiomics6 python=3.8
	conda activate multiomics6
	module load gcc/9.2.0
	pip install hic-straw
	conda install numpy
	conda install matplotlib
	conda install seaborn
	conda install requests
	conda install h5py
	conda install scipy
	pip install pyBigWig
	pip install fanc

#STEP2: GET THE HIC DATA
I have two shell files that make the file structure and generate the .mat files to work with
	sbatch makeFileStructure.sh
	sbatch getData.sh

This will produce the pearsons matrices, translated pearosns matrices, eigenvector, and hic observed matrices. It is running now, but is taking a long time computing the Pearsons matrices for the resolution 10,000 matrices
I am also running a version that computed the translated pearsons for all except the resolution 10,000

#STEP3: 
sbatch getStructuredData.sh

The Tensorlab Matlab package needs to be installed. Also, the tensorlab file "mtkrprod.m" needs to be edited before running low rank decomposition on high resolution data. Specficially, the line
largescale = mem > 2e9 || (isstruct(data) && isempty(data.matrix));
needs to be changed to 
largescale = mem > 2e15 || (isstruct(data) && isempty(data.matrix));

#STEP4: A/B Compartment Analysis
sbatch getOtherData.sh
sbatch getImages.sh

#STEP5: Insulation Score Analysis
getInsulationData.sh
sbatch getInsulationImages.sh


