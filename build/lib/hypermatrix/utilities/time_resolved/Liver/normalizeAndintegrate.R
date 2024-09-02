# Load required libraries
library(Seurat)
library(MASS)
library(DoubletFinder) 

# Function to add a column for the percent of mitochondrial genes ratio
add_mito_ratio_column <- function(seurat_obj, gene_pattern = "^MT-") {
  # Calculate percentage of mitochondrial genes if not already done
  if (!"percent.mt" %in% colnames(seurat_obj)) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = gene_pattern, col.name = "percent.mt")
  }
  
  return(seurat_obj)
}

# Function to print the max mitochondrial gene ratio for a sample
print_max_mito_ratio <- function(seurat_obj, sample_name) {
  # Extract the max mitochondrial gene ratio for the specified sample
  max_mito_ratio <- max(seurat_obj$percent.mt[seurat_obj$sample == sample_name], na.rm = TRUE)
  
  # Print the result
  cat("Max Mitochondrial Gene Ratio for", sample_name, ":", max_mito_ratio, "\n")
}

# Example usage:

# Assuming you have a Seurat object named 'seurat_obj'
# Add mitochondrial gene ratio column
seurat_obj_Het_O <- add_mito_ratio_column(seurat_obj_Het_O)

# Print the max mitochondrial gene ratio for a sample (replace 'YourSampleName' with the actual sample name)
print_max_mito_ratio(seurat_obj_Het_O, "YourSampleName")

print_max_mito_ratio <- function(seurat_obj) {
  # Identify unique sample identifiers
  sample_ids <- unique(seurat_obj$orig.ident)
  
  # Initialize variables
  max_ratio <- 0
  max_sample <- ""
  
  # Iterate through samples
  for (sample_id in sample_ids) {
    # Calculate max mitochondrial gene ratio for the current sample
    max_mito_ratio <- max(seurat_obj$percent.mt[seurat_obj$orig.ident == sample_id], na.rm = TRUE)
    
    # Check if the current sample has a higher ratio
    if (max_mito_ratio > max_ratio) {
      max_ratio <- max_mito_ratio
      max_sample <- sample_id
    }
  }
  
  # Print the result
  if (max_ratio > 0) {
    cat("Sample with Max Mitochondrial Gene Ratio:", max_sample, "\n")
    cat("Max Mitochondrial Gene Ratio:", max_ratio, "\n")
  } else {
    cat("No sample found.\n")
  }
}

# Example usage:

# Assuming you have a Seurat object named 'seurat_obj_Het_O'
# Add mitochondrial gene ratio column if not already done
if (!"percent.mt" %in% colnames(seurat_obj_Het_O)) {
  seurat_obj_Het_O <- PercentageFeatureSet(seurat_obj_Het_O, pattern = "^MT-", col.name = "percent.mt")
}

# Print the max mitochondrial gene ratio and corresponding sample
print_max_mito_ratio(seurat_obj_Het_O)

print_max_mito_ratio(seurat_obj_Het_O)

head(seurat_obj_Het_O)

# Assuming you have a Seurat object named 'seurat_obj'
# Add mitochondrial gene ratio column if not already done
if (!"percent.mt" %in% colnames(seurat_obj)) {
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
}

# Function to filter cells based on gene number and mitochondrial gene ratio
filter_cells <- function(seurat_obj, gene_threshold, mito_ratio_threshold) {
  # Identify cells to keep based on gene number
  keep_cells <- which(seurat_obj$nFeature_RNA > gene_threshold)
  
  # Calculate percentage of mitochondrial genes if not already done
  if (!"percent.mt" %in% colnames(seurat_obj)) {
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
  }
  
  # Further filter cells based on mitochondrial gene ratio
  keep_cells <- keep_cells[seurat_obj$percent.mt[keep_cells] <= mito_ratio_threshold]
  
  # Subset Seurat object to keep only the selected cells
  seurat_obj <- seurat_obj[, keep_cells]
  
  return(seurat_obj)
}

# Function to perform doublet detection
detect_doublets <- function(seurat_obj) {
  doublet_scores <- doubletFinder(seurat_obj)
  seurat_obj <- subset(seurat_obj, cells = which(doublet_scores$doublet == FALSE))
  
  return(seurat_obj)
}

# Load and preprocess individual datasets

path <- '/home/dwk681/workspace/CRA004660/Liver/Het-O/CRR403693_cellranger_output/outs/'
dataset_Het_O <- Read10X(paste(path, "filtered_feature_bc_matrix", sep = "/"))

path <- '/home/dwk681/workspace/CRA004660/Liver/Iso-Y/CRR403690_cellranger_output/outs/'
dataset_Iso_Y <- Read10X(paste(path, "filtered_feature_bc_matrix", sep = "/"))

path <- '/home/dwk681/workspace/CRA004660/Liver/Iso-O/CRR403691_cellranger_output/outs/'
dataset_Iso_O <- Read10X(paste(path, "filtered_feature_bc_matrix", sep = "/"))

path <- '/home/dwk681/workspace/CRA004660/Liver/Het-Y/CRR403692_cellranger_output/outs/'
dataset_Het_Y <- Read10X(paste(path, "filtered_feature_bc_matrix", sep = "/"))

# Create Seurat objects for each dataset
seurat_obj_Het_O <- CreateSeuratObject(counts = dataset_Het_O)
seurat_obj_Iso_Y <- CreateSeuratObject(counts = dataset_Iso_Y)
seurat_obj_Iso_O <- CreateSeuratObject(counts = dataset_Iso_O)
seurat_obj_Het_Y <- CreateSeuratObject(counts = dataset_Het_Y)

# Filter cells based on gene number and mitochondrial gene ratio
seurat_obj_Het_O <- filter_cells(seurat_obj_Het_O, 200, 0.05)
seurat_obj_Iso_Y <- filter_cells(seurat_obj_Iso_Y, 200, 0.05)
seurat_obj_Iso_O <- filter_cells(seurat_obj_Iso_O, 200, 0.05)
seurat_obj_Het_Y <- filter_cells(seurat_obj_Het_Y, 200, 0.05)

# Doublet detection
seurat_obj_Het_O <- detect_doublets(seurat_obj_Het_O)
seurat_obj_Iso_Y <- detect_doublets(seurat_obj_Iso_Y)
seurat_obj_Iso_O <- detect_doublets(seurat_obj_Iso_O)
seurat_obj_Het_Y <- detect_doublets(seurat_obj_Het_Y)


# Normalize each dataset using SCTransform
seurat_obj_Het_O <- SCTransform(seurat_obj_Het_O)
seurat_obj_Iso_Y <- SCTransform(seurat_obj_Iso_Y)
seurat_obj_Iso_O <- SCTransform(seurat_obj_Iso_O)
seurat_obj_Het_Y <- SCTransform(seurat_obj_Het_Y)

# Integrate the datasets
seurat_objs <- list(seurat_obj_Het_O, seurat_obj_Iso_Y, seurat_obj_Iso_O, seurat_obj_Het_Y)
anchors <- FindIntegrationAnchors(seurat_objs)
seurat_integrated <- PrepSCTIntegration(seurat_objs, anchors = anchors)

# Integrate the datasets using Canonical Correlation Analysis (CCA)
seurat_integrated <- IntegrateData(seurat_integrated)

# Extract the expression matrix, gene names, and cell names
expression_matrix <- GetAssayData(seurat_integrated)
gene_names <- rownames(expression_matrix)
cell_names <- colnames(expression_matrix)

# Save the expression matrix, gene names, and cell names as a MATLAB file
mat_file <- "expression_matrix.mat"
write.matrix(expression_matrix, mat_file, varnames = list(genes = gene_names, cells = cell_names))

# Perform dimensionality reduction (e.g., PCA)
seurat_integrated <- RunPCA(seurat_integrated)

# Cluster and visualize the integrated data
seurat_integrated <- FindNeighbors(seurat_integrated)
seurat_integrated <- FindClusters(seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated)

markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,  # Consider only positive markers
  min.pct = 0.1,    # Minimum percentage of cells expressing the gene
  thresh.use = 0.25,  # Threshold for the difference of the means in the natural log scale
  test.use = "wilcox",  # Wilcoxon rank-sum test
  logfc.threshold = 0.5,  # Log-fold change threshold
  min.diff.pct = 0.1,  # Minimum percentage of cells expressing a gene in one cluster compared to others
  random.seed = 123,  # Set a random seed for reproducibility
  return.thresh = 0.05  # Return only markers with adjusted p-value below this threshold
)

# Access the results
significant_markers <- markers[markers$adj.P.Val < 0.05 & abs(markers$logFC) > 0.5, ]

# Print or save the results
print(significant_markers)


