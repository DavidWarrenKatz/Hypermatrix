#!/bin/bash

# Clone the repository from Bitbucket
git clone https://bitbucket.org/dnaase/bisulfitehic.git

# Unzip the cloned repository (assumes the zip file exists in the current directory)
unzip dnaase-bisulfitehic-5dbb6ce7bc3c.zip

# Navigate into the unzipped directory
cd dnaase-bisulfitehic-5dbb6ce7bc3c

# Source the install.sh script
source ./install.sh

# Create a new conda environment named 'bisulfitehic' with Python 3.6
conda create -n bisulfitehic python=3.6 -y

# Activate the conda environment
source activate bisulfitehic

# Install the required Python packages
pip install numpy
pip install pysam

# Set the JAVA_HOME environment variable
export JAVA_HOME=$(/usr/libexec/java_home)

echo "Setup complete. The 'bisulfitehic' environment is activated."

