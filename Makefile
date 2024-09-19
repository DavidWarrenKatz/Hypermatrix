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
CONDA_ENV_PATH = $(if $(CUSTOM_ENV_PATH),$(CUSTOM_ENV_PATH),$(CONDA_BASE)/envs)

.PHONY: all env install clean backup test

all: env install

env:
	@echo "Checking if Conda environment '$(ENV_NAME)' exists..."
	@if conda env list | grep -q "^$(ENV_NAME)\s"; then \
		echo "Environment '$(ENV_NAME)' already exists at $(CONDA_ENV_PATH)/$(ENV_NAME)"; \
	else \
		echo "Creating Conda environment '$(ENV_NAME)'..."; \
		conda env create -f $(YML_FILE) -n $(ENV_NAME); \
	fi
	@echo "Ensuring 'r-base' and 'gcc' are installed in the environment..."
	$(CONDA_ACTIVATE) && conda install -y r-base gcc

install: env
	@echo "Activating environment and installing Python dependencies..."
	$(CONDA_ACTIVATE) && \
	pip install . && \
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
	@echo "Removing Conda environment '$(ENV_NAME)'..."
	conda remove --name $(ENV_NAME) --all --yes

backup:
	@echo "Backing up Conda environment to '$(ENV_NAME)_backup.yml'..."
	$(CONDA_ACTIVATE) && conda env export > $(ENV_NAME)_backup.yml

test:
	@echo "Testing Conda environment '$(ENV_NAME)'..."
	$(CONDA_ACTIVATE) && \
	python -c "import sys; print('Python version:', sys.version)" && \
	R --version && \
	gcc --version
