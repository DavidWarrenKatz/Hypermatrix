# Hypermatrix

Advancements in single-cell multi-omics technologies have enabled the simultaneous measurement of various omics modalities within individual cells. Integrating multi-omics data while preserving the interaction information between different modalities remains an open challenge. Traditional methods lose critical interaction information by applying matrix methods. To address this, this research project proposes a Non-Negative Tensor Factorization (NTF) model for multi-omics integration called HYPERMATRIX.

This software derives cell-type factors as well as A/B compartment factors from integrated bulk epigenetic data and integrated single-cell Hi-C and methylation data.

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

4. Navigate to the utilities directory:

    ```bash
    cd hypermatrix/utilities
    ```

5. Make the `getData.sh` script executable:

    ```bash
    chmod +x getData.sh
    ```

6. Run the `getData.sh` script to get the necessary data:

    ```bash
    ./getData.sh
    ```

## Usage

After setting up the environment and downloading the necessary data, you can start using the HYPERMATRIX tool. Provide detailed usage instructions here, including any example commands or scripts that users can run.

```bash
# Example command
python run_hypermatrix.py --input data/input_file --output results/output_file
