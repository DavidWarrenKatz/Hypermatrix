# To-Do List

- [ ] Task 1: ABcluster
      incorporate KR normalization and ICE normalization functions into single-cell pipeline.
      make command command-line callable
      make function to find differentially expressed regions to display for KR normalized runs.
      add the clustering cells pipeline code.
- [ ] Task 2: differentiate_chromosomes
      process the iterchromosomal contacts seperately with a flag so tha thte pipeline does not have to be run on the genome-wide tensor. 
- [ ] Task 3: Noise Reduction before loop calling with Hiccups
- [ ] Task 4: Finish Bulk Pipeline
- [ ] Task 5: Finish synthetic pipeline and upload.
      create file to test cumulant pipeline on synthetic data
      create file to demonstrate how there are in fact two compartments, and the Hi-C data cannot be explain by arbitray even numbered compartments. 
- [ ] Task 6: Make everything run in parallel where possible
- [ ] Task 7: Compare Hypermatrix with scHiCLuster
- [ ] Task 8: Compare Hypermatrix with Higashi
- [ ] Task 9: Do through Dixon et al again and look at their omics data
- [ ] Task 10: Do some analysis on genomic factors
- [ ] Guided factorization
        Divide the genome into bins. Make the following category vectors based on data.
Isochore - Classify bins into categories L1, L2, H1, H2, and H3 based on GC contents - <38%, 38%-42%, 42%-47%, 47%-52%, and >52%, respectively
Replication timing - Classify bins into categories E0, E1, L1, L2, and L3 based on replication timing.
Make a matrix with entry (i,j) equal to 1 with i and j are in the same category, and 0 if they are in different categories.
- [ ] 

      
## Annotated Bibliography
- [ ] Yao, Z., Van Velthoven, C. T., Nguyen, T. N., Goldy, J., Sedeno-Cortes, A. E., Baftizadeh, F., ... & Zeng, H. (2021). A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell, 184(12), 3222-3241.
- [ ] Zeisel, A., Muñoz-Manchado, A. B., Codeluppi, S., Lönnerberg, P., La Manno, G., Juréus, A., ... & Linnarsson, S. (2015). Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science, 347(6226), 1138-1142.
- [ ] Luo, C., Keown, C. L., Kurihara, L., Zhou, J., He, Y., Li, J., ... & Ecker, J. R. (2017). Single-cell methylomes identify neuronal subtypes and regulatory elements in mammalian cortex. Science, 357(6351), 600-604.
- [ ] Luo, C., Liu, H., Xie, F., Armand, E. J., Siletti, K., Bakken, T. E., ... & Ecker, J. R. (2022). Single nucleus multi-omics identifies human cortical cell regulatory genome diversity. Cell genomics, 2(3).
- [ ] Luo, C., Rivkin, A., Zhou, J., Sandoval, J. P., Kurihara, L., Lucero, J., ... & Ecker, J. R. (2018). Robust single-cell DNA methylome profiling with snmC-seq2. Nature communications, 9(1), 3824.
- [ ] Liu, H., Zeng, Q., Zhou, J., Bartlett, A., Wang, B. A., Berube, P., ... & Ecker, J. R. (2023). Single-cell DNA methylome and 3D multi-omic atlas of the adult mouse brain. Nature, 624(7991), 366-377.
- [ ] Lee, D. S., Luo, C., Zhou, J., Chandran, S., Rivkin, A., Bartlett, A., ... & Ecker, J. R. (2019). Simultaneous profiling of 3D genome structure and DNA methylation in single human cells. Nature methods, 16(10), 999-1006.
- [ ] Zhang, R., Zhou, T., & Ma, J. (2022). Ultrafast and interpretable single-cell 3D genome analysis with Fast-Higashi. Cell systems, 13(10), 798-807.
- [ ] Tan, L., Xing, D., Chang, C. H., Li, H., & Xie, X. S. (2018). Three-dimensional genome structures of single diploid human cells. Science, 361(6405), 924-928.
      method for calling single-cell compartments and seperating homologous chromosomes. sc  compartment value of bin i is the weighted average CpG density in the reference genome of all the bins that bin i contacts in the single cell contact map. 
- [ ] Nagano, T., Lubling, Y., Várnai, C., Dudley, C., Leung, W., Baran, Y., ... & Tanay, A. (2017). Cell-cycle dynamics of chromosomal organization at single-cell resolution. Nature, 547(7661), 61-67.
      sc  compartment value of bin i is the weighted average A/B compartment value of bulk Hi-C of all the bins that bin i contacts in the single cell contact map. 
- [ ] Zhou, J., Ma, J., Chen, Y., Cheng, C., Bao, B., Peng, J., ... & Ecker, J. R. (2019). Robust single-cell Hi-C clustering by convolution-and random-walk–based imputation. Proceedings of the National Academy of Sciences, 116(28), 14011-14018.
      sc  compartment value of bin i  is the dot product of imputed Hi-C contact maps with CpG density vector.
- [ ] Xie, W. J., Meng, L., Liu, S., Zhang, L., Cai, X., & Gao, Y. Q. (2017). Structural modeling of chromatin integrates genome features and reveals chromosome folding principle. Scientific reports, 7(1), 2818.
      The authors show that average CpG frequency from reference genome can serve as a good proxy of A/B compartments in bulk Hi-C data. The large majority of partially methylated domains as called by scanning the reference genome are in the B compartment as called by principal component analysis. 
- [ ] Benjamini, Y., & Speed, T. P. (2012). Summarizing and correcting the GC content bias in high-throughput sequencing. Nucleic acids research, 40(10), e72-e72.
- [ ] Welch, J. D., Kozareva, V., Ferreira, A., Vanderburg, C., Martin, C., & Macosko, E. Z. (2019). Single-cell multi-omic integration compares and contrasts features of brain cell identity. Cell, 177(7), 1873-1887.
      3-fold tensor method
- [ ] Stevens, T. J., Lando, D., Basu, S., Atkinson, L. P., Cao, Y., Lee, S. F., ... & Laue, E. D. (2017). 3D structures of individual mammalian genomes studied by single-cell Hi-C. Nature, 544(7648), 59-64.
- [ ] Kind, J., Pagie, L., Ortabozkoyun, H., Boyle, S., De Vries, S. S., Janssen, H., ... & van Steensel, B. (2013). Single-cell dynamics of genome-nuclear lamina interactions. Cell, 153(1), 178-192.
