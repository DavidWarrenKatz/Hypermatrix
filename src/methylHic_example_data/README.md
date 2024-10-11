# Single Cell Data

The `SRR_Acc_List_sc.txt` contains the SRAs of single cells. There are **59 SRA files** that contain the 4-multiplexed reads from **150 cells**.

### Primed Cells

- **Number of cells**: 103
- **SRR IDs**: SRR7770822 - SRR7770868
- **Culture Condition**: Serum + LIF

### Naive Cells

- **Number of cells**: 47
- **SRR IDs**: SRR7770869 - SRR7770880
- **Culture Condition**: 2i + LIF (except for Naive 1-4 which is serum + LIF, possibly a mistake in SRA SRR7770869)


The script download_and_
extracts 118 fastq files, two for each of the 59 SRA files. 
after demultiplexing, we ar eleft with 150 bam files

---

# Bulk Data

The `SRR_Acc_List_bulk.txt` contains the SRAs of the bulk data.
 
