# Makefile for setting up the Hypermatrix environment

.ONESHELL:
SHELL = /bin/bash
ENV_NAME = hypermatrix
YML_FILE = hypermatrix.yml
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate $(ENV_NAME)

.PHONY: all env install clean

all: env install

# Create Conda environment
env:
	@echo "Creating Conda environment..."
	conda env create -f $(YML_FILE) --yes

# Install Python dependencies using pip
install:
	@echo "Activating environment and installing Python dependencies..."
	$(CONDA_ACTIVATE) && pip install .

# Clean up the environment
clean:
	@echo "Removing Conda environment..."
	conda remove --name $(ENV_NAME) --all --yes
