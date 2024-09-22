# Hypermatrix

Hypermatrix is a command-line tool designed for the integration of multi-omics data, as well as general epigenetic data analysis using tensor techniques. Advancements in single-cell multi-omics technologies have enabled the simultaneous measurement of various omics modalities within individual cells. Integrating multi-omics data while preserving the interaction information between different modalities remains an open challenge. Traditional matrix methods lose critical interaction information. To address this, the 'Hypermatrix' software implements a Non-Negative Tensor Factorization (NTF) pipeline for multi-omics integration.

# Key Commands
**ABcluster**: The 'ABcluster' command processes single-cell CpG methylation and single-cell chromosome conformation data to perform cell-type clustering and single-cell A/B compartment identification. This command is compatible with data from techniques like sn-m3C-seq (described by Dong-Sung Lee et al. in Nature Methods, 2019), scMethyl-HiC (described by Guoqiang Li et al. in Nature Methods, 2019) and the single-cell version of NOMe-HiC (described by Hailu Fu et al. in Genome Biology, 2023). 'ABcluster' is particularly useful for tracking the diversity of A/B compartments within cells of the same type. This tool can be used to test the hypothesis that the diversity in A/B compartments within a cell-type increases with the age of an organism.

**TADcluster** *(in progress)*: The 'TADcluster' command processes single-cell CpG methylation and single-cell chromosome conformation data to perform cell-type clustering and single-cell TAD boundary detection. 

**differentiate_chromosomes** *(in progress)*: The 'differentiate_chromosomes' command processes Hi-C data, optionally combined with other epigenetic modalities, to produce distinct A/B compartment calls for each homologous chromosome. Unlike 'ABcluster', which relies solely on intrachromosomal contacts, 'differentiate_chromosomes' uses both intra- and interchromosomal contacts. The folding of chromatin in 3D space is an important part of gene regulation, ensuring that certain genes are transcribed simultaneously and that their transcripts are spatially close for further processing. The 'differentiate_chromosomes' command can be used to test the hypothesis that diploid cells use two distinct regulatory configurations for each homologous chromosome, with each 3D conformation likely being mutually exclusive to the other. For example, for one folding program to position certain genes in an active transcriptional hub, it may have to seperated other genes that also need to be clustered. Diploidy provides a solution to this problem. Diploidy offers more than just a backup copy of each chromosome; it provides cells with additional regulatory complexity by enabling the use of two mutually exclusive folding programs simultaneously. Additionally, the 'differentiate_chromosomes' command measures the degree to which each chromosome is associated with the nuclear lamina.

## Visualization

Below is the heatmap comparing the bulk eigenvectors of GM12878 and IMR90 with single-cell compartment calls. The first 38 cells are GM12878.

<div style="text-align: center;">
  <table style="margin: 0 auto;">
    <tr>
      <td><img src="output_files/AB_compartment_heatmap_ch10_example.png" alt="Figure 1" width="500"></td>
      <td><img src="output_files/AB_compartment_heatmap_ch4_example.png" alt="Figure 2" width="550"></td>
    </tr>
  </table>
</div>

## Installation

### Prerequisites

Ensure you have `conda` installed. If not, you can install it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). The packages `pyBigWig`, `h5py`, `hic-straw`, and `scHiCluster` are required and will be installed during the process below. If any of the required packages are not properly installed, you may need to install them manually.

### Steps

1. Clone the repository:

    ```bash
    git clone https://github.com/DavidWarrenKatz/hypermatrix.git
    ```

2. Navigate into the cloned directory:

    ```bash
    cd hypermatrix
    ```

3. Install the Hypermatrix tool and its dependencies:

    ```bash
    make install
    ```

4. Activate the Conda environment:

    ```bash
    conda activate hypermatrix
    ```

5. Verify the installation:

    ```bash
    hypermatrix --version
    ```

    This command should display something like:
    
    ```
    Hypermatrix version 0.1 - A tool for integrating multi-omics data and epigenetic analysis using advanced tensor techniques.
    ```

### ABcluster Command

The `ABcluster` command is used to perform single-cell A/B compartment analysis and identify cell-type clusters by integrating single-cell CpG methylation data and Hi-C data. This command uses one or both data modalities depending on the user's input.

#### General Syntax

```bash
hypermatrix ABcluster --methy <path_to_methylation_directory> --hic <path_to_hic_directory> --output_dir <output_directory>
```

#### Input Parameters

- **`--methy <path_to_methylation_directory>`**: This specifies the directory containing the single-cell CpG methylation files. These files must be named following the pattern `<prefix>.bw`, where `<prefix>` is a unique identifier for each sample.
  
- **`--hic <path_to_hic_directory>`**: This specifies the directory containing the single-cell Hi-C files. These files must follow the naming pattern `<prefix>.hic`, where `<prefix>` matches the one used in the methylation files for proper integration.

- **`--output_dir <output_directory>`**: This specifies the directory where the output results will be stored. The output directory will contain a file `cell_type_clusters.txt` for each cell-type, and a file `prefix_ab_compartments.txt` for the A/B compartment call for each cell.

#### Configurable Parameters

The parameters for the `ABcluster` command are listed in the file `hypermatrix/config.py`, where they can be adjusted.

#### Usage Recommendations

It is recommended to use both the `--methy` and `--hic` flags together when both types of data are available. 
  
- **Using only the `--methy` flag**: If only the methylation data is available, the software will generate A/B compartment calls and cell-type clusters based solely on the single-cell CpG methylation data.

- **Using only the `--hic` flag**: Similarly, if only Hi-C data is available, the software will generate A/B compartment calls and cell-type clusters based solely on the single-cell Hi-C data. The results obtained from using only Hi-C data can be compared to those generated by scHiCluster and FastHigashi software.


### **Preprocess Command**

The `preprocess` command processes BAM files to prepare them for use in commands like `ABcluster`. If run with the `--nomehic` flag, the BAM files are assumed to be from the scNOMe-HiC technique and are processed accordingly. This will produce Hi-C and methylation files in the correct format for the `ABcluster` command.

#### General Syntax

```bash
hypermatrix preprocess --nomehic --input_dir <path_to_bam_directory> --output_dir <path_to_output_directory> --ref_genome <path_to_reference_genome>
```

#### Input Parameters

- **`--input_dir <path_to_bam_directory>`**: Specifies the directory containing the indexed BAM files for processing. It is assumed that all BAM files are indexed (i.e., corresponding `.bai` files are present).

- **`--output_dir <path_to_output_directory>`**: Specifies the output directory. A subfolder for Hi-C, methylation, and accessibility will be created.
  
- **`--ref_genome <path_to_reference_genome>`**: Specifies the reference genome. Only options right now are hg19 and hg38. 

## Example Usage of Hypermatrix Software

First, download single-cell Methyl-HiC fastq files from SRA under accession number SRP159191 using the `prefetch` and `fastq-dump` utilities from the SRA toolkit. An example scripts for doing this are located in the directory hypermatrix/src/methylHic_example_data. Thne, run the following command to preprocess the Methyl-HiC data. After preprocessing, run `ABcluster` cluster command.

```bash
hypermatrix preprocess --nomehic --input_dir /path/to/fastq_folder ---output_dir /path/to/output -ref hg19
hypermatrix ABcluster --methy /path/to/output/methy --hic /path/to/output/hic --output_dir <output_directory> --res 1000000
```

## Contact

For any questions or issues, please contact davidkatz02@gmail.com.
