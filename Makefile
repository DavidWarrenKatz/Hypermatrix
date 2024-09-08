.ONESHELL:
SHELL = /bin/bash
ENV_NAME = hypermatrix
YML_FILE = hypermatrix.yml
CONDA_BASE = $(shell conda info --base)
CONDA_ACTIVATE = source $(CONDA_BASE)/etc/profile.d/conda.sh && conda activate $(ENV_NAME)

.PHONY: all env install clean

all: env install

env:
	@echo "Creating Conda environment..."
	conda env create -f $(YML_FILE)

install:
	@echo "Activating environment and installing Python dependencies..."
	$(CONDA_ACTIVATE) && pip install . && pip install --use-pep517 fanc hic-straw

clean:
	@echo "Removing Conda environment..."
	conda remove --name $(ENV_NAME) --all --yes
