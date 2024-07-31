# Hypermatrix

Hypermatrix is a command-line tool designed for the integration of multi-omics data. The 'ABcluster' command inputs single-cell CpG methylation and chromosome conformation information to perform cell-type clustering, A/B compartment calls, and TAD boundary calls for each cell. Example data to input into this command is Dong-Sung Lee et al. (Nature Methods, 2019) and Hailu Fu et al. (Genome Biology, 2023). The 'differentiate_chromosomes' command differentes between the homologous chromosomes and determines if B compartments are lamina-associated.  

Advancements in single-cell multi-omics technologies have enabled the simultaneous measurement of various omics modalities within individual cells. Integrating multi-omics data while preserving the interaction information between different modalities remains an open challenge. Traditional methods lose critical interaction information by applying matrix methods. To address this, this research project proposes a Non-Negative Tensor Factorization (NTF) model for multi-omics integration called HYPERMATRIX.

## Installation

### Prerequisites

Ensure you have `conda` installed. If not, you can install it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Steps

1. Clone the repository:

    ```bash
    git clone https://github.com/DavidWarrenKatz/hypermatrix.git
    ```

2. Create the conda environment:

    ```bash
    conda env create -f hypermatrix/bin/hypermatrix.yml
    ```

3. Activate the environment:

    ```bash
    conda activate hypermatrix
    ```

4. As a first demonstration, derives A/B compartments integrated bulk epigenetic data and and synthetic data to illustrate the non-negative tensor decomposition method.
 Navigate to the utilities directory:

    ```bash
    cd hypermatrix/utilities
    ```

5. Make the `runPipeline_bulkdata.sh` script executable:

    ```bash
    chmod +x runPipeline_bulkdata.sh
    ```

6. Run the `runPipeline_bulkdata.sh` script to get the necessary data:

    ```bash
    ./runPipeline_bulkdata.sh
    ```

## Usage

After setting up the environment and downloading the necessary data, you can start using the HYPERMATRIX tool. Provide detailed usage instructions here, including any example commands or scripts that users can run.

```bash
# Example command
python run_hypermatrix.py --input data/input_file --output results/output_file
