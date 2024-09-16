.ONESHELL:
SHELL = /bin/bash
ENV_NAME = hypermatrix
YML_FILE = hypermatrix/installation/hypermatrix.yml
CONDA_BASE = $(shell conda info --base)
CONDA_ACTIVATE = source $(CONDA_BASE)/etc/profile.d/conda.sh && conda activate $(ENV_NAME)

# Uncomment the following line to specify a custom environment directory (e.g., /path/to/custom_envs)
# CUSTOM_ENV_PATH = /home/dkatz/conda_envs
#
# Use default environment path unless a custom one is provided
# Include r-base 
# Include gcc (this is an important requirement)
# Add env test to backup installations	
CONDA_ENV_PATH = $(if $(CUSTOM_ENV_PATH),$(CUSTOM_ENV_PATH),$(CONDA_BASE)/envs)

.PHONY: all env install clean

all: env install

env:
	@echo "Checking if Conda environment exists..."
	@if conda env list | grep -q "^$(ENV_NAME)\s"; then \
		echo "Environment '$(ENV_NAME)' already exists at $(CONDA_ENV_PATH)/$(ENV_NAME)"; \
	else \
		echo "Creating Conda environment '$(ENV_NAME)'..."; \
		conda env create -f $(YML_FILE) -n $(ENV_NAME); \
	fi

install: env
	@echo "Activating environment and installing Python dependencies..."
	$(CONDA_ACTIVATE) && pip install . && \
	pip install --use-pep517 fanc hic-straw && \
	pip install git+https://github.com/zhoujt1994/scHiCluster.git || \
	( \
		echo "Initial install failed, falling back to loading modules and setting flags..." && \
		module load gcc/9.3.0 R && \
		echo "Activating environment and retrying installation..." && \
		$(CONDA_ACTIVATE) && \
		export CFLAGS="-std=c99" && \
		pip install . && \
		pip install --use-pep517 fanc hic-straw && \
		pip install git+https://github.com/zhoujt1994/scHiCluster.git \
	)

shell:
	@echo "Launching shell with activated environment..."
	$(CONDA_ACTIVATE) && /bin/bash

clean:
	@echo "Removing Conda environment..."
	conda remove --name $(ENV_NAME) --all --yes






