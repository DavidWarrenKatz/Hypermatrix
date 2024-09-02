#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=1-0:00:00
#SBATCH --ntasks=5
#SBATCH --mem=150000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetNMFPlot.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetNMFPlot.out


python << EOF
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import NMF
from tensorly.decomposition import non_negative_parafac
import tensorly as tl
import anndata

file_path = '/home/dwk681/workspace/CRA004660/Liver/integrated_and_normalized_annData.h5ad'

try:
    integrated_adata = anndata.read_h5ad(file_path)
    print(f"AnnData object successfully loaded from: {file_path}")
except Exception as e:
    print(f"Error: {e}")

data_matrix = integrated_adata.X.toarray()



def compute_nmf_tsne_and_save_plot(matrix, rank, title, save_path, n_components=2, perplexity=10, random_state=None):
    # Compute NMF for the given rank
    nmf_model = NMF(n_components=rank, init='random', random_state=42, max_iter=200)
    nmf_components = nmf_model.fit_transform(matrix)

    # Assuming you want to visualize the first mode of NMF components
    ntf_components = nmf_components[:, 0]

    # Perform t-SNE on NMF components
    tsne_result = TSNE(n_components=n_components, perplexity=perplexity, random_state=random_state).fit_transform(ntf_components.reshape(-1, 1))
    # Plot the t-SNE result
    fig, ax = plt.subplots(figsize=(8, 8))
    scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1])
    ax.set_title(title, fontsize=30)
    ax.set_xlabel('t-SNE Component 1', fontsize=20)
    ax.set_ylabel('t-SNE Component 2', fontsize=20)

    # Save the figure
    plt.savefig(save_path)

    return scatter




def make_matrix_non_negative(matrix):
    min_val = np.min(matrix)
    if min_val < 0:
        translated_matrix = matrix - min_val
        return translated_matrix
    else:
        return matrix

save_path = '/home/dwk681/workspace/CRA004660/Liver/tsneFigure.png'
data_matrix = make_matrix_non_negative(data_matrix)  
compute_nmf_tsne_and_save_plot(data_matrix, rank=5, title='tSNE of NMF Components', save_path=save_path, n_components=2, perplexity=5, random_state=42)



EOF
