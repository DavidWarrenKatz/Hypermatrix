# hypermatrix
 Advancements in single-cell multi-omics technologies have enabled the simultaneous measurement of various omics modalities within individual cells. Integrating multi-omics data while preserving the interaction information between different modalities remains an open challenge. Traditional methods lose critical interaction information by applying matrix methods. To address this, this research project proposes a Non-Negative Tensor Factorization (NTF) model for multi-omics integration called HYPERMATRIX. 

This preliminary version derived cell-type factors as well as A/B compartment facotrs from integrated bulk epigenetic data and integrated single-cell Hi-C and methylation data. 

Follow the following steps:
git clone https://github.com/DavidWarrenKatz/hypermatrix.git

cd hypermatrix/utilities

chmod +x getData.sh 

./getData.sh

This file downloads bulk Hi-C files from a reference Hi-C file. All reference files are GM12878 by defauls, but can be easily changed as needed.
