#!/bin/bash
#SBATCH --job-name=bwa_indexing
#SBATCH --output=bwa_indexing_%j.log
#SBATCH --error=bwa_indexing_%j.err
#SBATCH --time=12:00:00    # Set the max time to 12 hours
#SBATCH --mem=32G          # Request 32GB of memory (adjust as necessary)
#SBATCH --ntasks=1         # Request a single task
#SBATCH --cpus-per-task=8  # Use 8 CPUs for faster indexing (adjust as necessary)

#./download_example_data.sh
./setup_bhmem.sh
