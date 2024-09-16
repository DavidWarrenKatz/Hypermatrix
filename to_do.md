# To-Do List

## Task 1: ABcluster
- [ x ] include utilities/**/* in setup.py 
- [ ] Incorporate KR normalization and ICE normalization functions into the single-cell pipeline.
- [ ] Make command-line callable.
- [ ] Create a function to find differentially expressed regions for KR normalized runs.
- [ ] Add the clustering cells pipeline code.
- [ ] Write the SNPS code

## Task 2: differentiate_chromosomes
- [ ] Process the interchromosomal contacts separately with a flag, so the pipeline does not have to run on the genome-wide tensor.

## Task 3: Noise Reduction
- [ ] Implement noise reduction before loop calling with Hiccups.

## Task 4: Finish Bulk Pipeline
- [ ] Complete the bulk pipeline.

## Task 5: Finish Synthetic Pipeline and Upload
- [ ] Create a file to test the cumulant pipeline on synthetic data.
- [ ] Create a file to demonstrate how there are, in fact, two compartments, and that the Hi-C data cannot be explained by arbitrary even-numbered compartments.

## Task 6: Parallelization
- [ ] Make everything run in parallel where possible.

## Task 7: Comparison with scHiCluster
- [ ] Compare Hypermatrix with scHiCluster.

## Task 8: Comparison with Higashi
- [ ] Compare Hypermatrix with Higashi.

## Task 9: Review Dixon et al.
- [ ] Go through Dixon et al. again and review their omics data.

## Task 10: Genomic Factors Analysis
- [ ] Perform analysis on genomic factors.

### Guided Factorization
- [ ] Divide the genome into bins. Create the following category vectors based on data:
  - **Isochore:** Classify bins into categories L1, L2, H1, H2, and H3 based on GC content:
    - L1: <38%
    - L2: 38%-42%
    - H1: 42%-47%
    - H2: 47%-52%
    - H3: >52%
  - **Replication Timing:** Classify bins into categories E0, E1, L1, L2, and L3 based on replication timing.
- [ ] Create a matrix with entry (i,j) equal to 1 if i and j are in the same category, and 0 if they are in different categories.

---

## Annotated Bibliography

- **Yao, Z., Van Velthoven, C. T., Nguyen, T. N., et al. (2021).**  
  A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. *Cell, 184(12), 3222-3241.*

- **Zeisel, A., Muñoz-Manchado, A. B., Codeluppi, S., et al. (2015).**  
  Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. *Science, 347(6226), 1138-1142.*

- **Luo, C., Keown, C. L., Kurihara, L., et al. (2017).**  
  Single-cell methylomes identify neuronal subtypes and regulatory elements in the mammalian cortex. *Science, 357(6351), 600-604.*

- **Luo, C., Liu, H., Xie, F., et al. (2022).**  
  Single nucleus multi-omics identifies human cortical cell regulatory genome diversity. *Cell Genomics, 2(3).*

- **Luo, C., Rivkin, A., Zhou, J., et al. (2018).**  
  Robust single-cell DNA methylome profiling with snmC-seq2. *Nature Communications, 9(1), 3824.*

- **Liu, H., Zeng, Q., Zhou, J., et al. (2023).**  
  Single-cell DNA methylome and 3D multi-omic atlas of the adult mouse brain. *Nature, 624(7991), 366-377.*

- **Lee, D. S., Luo, C., Zhou, J., et al. (2019).**  
  Simultaneous profiling of 3D genome structure and DNA methylation in single human cells. *Nature Methods, 16(10), 999-1006.*

- **Zhang, R., Zhou, T., & Ma, J. (2022).**  
  Ultrafast and interpretable single-cell 3D genome analysis with Fast-Higashi. *Cell Systems, 13(10), 798-807.*

- **Tan, L., Xing, D., Chang, C. H., et al. (2018).**  
  Three-dimensional genome structures of single diploid human cells. *Science, 361(6405), 924-928.*
  - Method for calling single-cell compartments and separating homologous chromosomes.
  - Single-cell compartment value of bin i is the weighted average CpG density in the reference genome of all the bins that bin i contacts in the single cell contact map.

- **Nagano, T., Lubling, Y., Várnai, C., et al. (2017).**  
  Cell-cycle dynamics of chromosomal organization at single-cell resolution. *Nature, 547(7661), 61-67.*
  - Single-cell compartment value of bin i is the weighted average A/B compartment value of bulk Hi-C of all the bins that bin i contacts in the single cell contact map.

- **Zhou, J., Ma, J., Chen, Y., et al. (2019).**  
  Robust single-cell Hi-C clustering by convolution-and random-walk–based imputation. *Proceedings of the National Academy of Sciences, 116(28), 14011-14018.*
  - Single-cell compartment value of bin i is the dot product of imputed Hi-C contact maps with CpG density vector.

- **Xie, W. J., Meng, L., Liu, S., et al. (2017).**  
  Structural modeling of chromatin integrates genome features and reveals chromosome folding principle. *Scientific Reports, 7(1), 2818.*
  - The authors show that average CpG frequency from the reference genome can serve as a good proxy for A/B compartments in bulk Hi-C data. The large majority of partially methylated domains as called by scanning the reference genome are in the B compartment as called by principal component analysis.

- **Benjamini, Y., & Speed, T. P. (2012).**  
  Summarizing and correcting the GC content bias in high-throughput sequencing. *Nucleic Acids Research, 40(10), e72-e72.*

- **Welch, J. D., Kozareva, V., Ferreira, A., et al. (2019).**  
  Single-cell multi-omic integration compares and contrasts features of brain cell identity. *Cell, 177(7), 1873-1887.*
  - 3-fold tensor method.

- **Stevens, T. J., Lando, D., Basu, S., et al. (2017).**  
  3D structures of individual mammalian genomes studied by single-cell Hi-C. *Nature, 544(7648), 59-64.*

- **Kind, J., Pagie, L., Ortabozkoyun, H., et al. (2013).**  
  Single-cell dynamics of genome-nuclear lamina interactions. *Cell, 153(1), 178-192.*





  Incorporate the other commands in readme

  ### Bulk Command

1. As a first demonstration, derive A/B compartments from integrated bulk epigenetic data and synthetic data to illustrate the non-negative tensor decomposition method. Navigate to the utilities directory:

    ```bash
    cd hypermatrix/utilities
    ```

2. Make the `runPipeline_bulkdata.sh` script executable:

    ```bash
    chmod +x runPipeline_bulkdata.sh
    ```

3. Run the `runPipeline_bulkdata.sh` script to get the necessary data:

    ```bash
    ./runPipeline_bulkdata.sh
    ```

### differentiate_chromosomes Command

To differentiate between the homologous chromosomes and determine if B compartments are lamina-associated:

```bash
python differentiate_chromosomes.py --input <path_to_input_file> --output <output_directory>
```

#### Arguments

- `--input`: Path to the input file containing chromosome data.
- `--output`: Directory where the output files will be saved.

#### Output

The output directory will contain:

- `homologous_chromosomes.txt`: Differentiated homologous chromosomes.
- `lamina_associated_B_compartments.txt`: B compartments that are lamina-associated.

#### Example Command

```bash
python differentiate_chromosomes.py --input data/chromosome_data.csv --output results/
```
