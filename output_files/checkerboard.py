import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def create_contact_matrix(labels):
    n = len(labels)
    contact_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if labels[i] == labels[j]:
                contact_matrix[i, j] = 1
    return contact_matrix

def save_heatmap_image(matrix, heatmap_filename):
    plt.imshow(matrix, cmap='hot', aspect='auto')
    plt.colorbar()
    plt.xlabel('genomic region')
    plt.ylabel('genomic region')
    plt.title('Contact Matrix Heatmap')
    plt.savefig(heatmap_filename)
    plt.close()

def save_scatterplot_image(labels, scatterplot_filename):
    n = len(labels)
    x = np.arange(n)
    y = np.zeros(n)
    plt.scatter(x, y, c=labels, cmap='Set1')
    plt.xlabel('genomic region')
    plt.ylabel('Cluster Label')
    plt.title('Cluster Labels Scatter Plot')
    plt.savefig(scatterplot_filename)
    plt.close()


def save_combined_image(matrix, labels, filename):
    # Create a figure with larger size
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    
    cmap = mcolors.ListedColormap(['white', 'black'])
    bounds = [0, 0.5, 1]  # Set color boundaries (0 for white, 1 for blue)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Plot heatmap with a larger font size for labels
    im = axes[0].imshow(matrix, cmap=cmap, aspect='auto', norm=norm)
    axes[0].set_xlabel('Genomic Region', fontsize=15)
    axes[0].set_ylabel('Genomic Region', fontsize=15)
    axes[0].set_title('Contact Matrix Heatmap', fontsize=20)

    # Customize the colorbar for the heatmap
    cbar = fig.colorbar(im, ax=axes[0], ticks=[0, 1])  # Set colorbar ticks to 0 and 1
    cbar.ax.set_yticklabels(['0', '1'])  # Set colorbar tick labels to '0' and '1'
    cbar.ax.tick_params(labelsize=10)
    

    # Plot scatter plot with larger dots and labels
    n = len(labels)
    x = np.arange(n)
    y = np.zeros(n)
    scatter = axes[1].scatter(x, y, c=labels, cmap='Set1', s=500)  # Increase s to make dots bigger
    axes[1].set_xlabel('Genomic Region', fontsize=15)
    #axes[1].set_ylabel('Cluster Label', fontsize=15)
    axes[1].set_title('Chromosome Colored By Cluster', fontsize=20)
    axes[1].tick_params(axis='y', labelleft=False)    

    legend_labels = ['Cluster 1', 'Cluster 2']

    handles, labels = scatter.legend_elements(prop='colors', alpha=0.9)
    legend = axes[1].legend(handles, legend_labels, loc='upper left')
    axes[1].add_artist(legend)

    # Customize the legend's font size and other properties
    legend.get_title().set_fontsize('20')  # Set the legend title font size
    for item in legend.get_texts():
        item.set_fontsize('15')  # Set the legend item (label) font size

    textbox_text = 'Each loci is assumed to interact \n' \
                   'equally with every other loci \n' \
                   'of its cluster and to not interact \n' \
                   'with any loci of the other cluster. \n' \
                   'This produces a perfect checkerboard.'
   
    axes[1].text(0.1, 0.1, textbox_text, fontsize=15, transform=axes[1].transAxes,
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

    # Customize the colorbar for the scatter plot
    #cbar = fig.colorbar(scatter, ax=axes[1])
    #cbar.ax.tick_params(labelsize=10)

    # Adjust spacing between subplots
    fig.tight_layout()

    # Save combined image with higher resolution (DPI)
    plt.savefig(filename, dpi=300)
    plt.close()

'''
def save_combined_image(matrix, labels, filename):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Plot heatmap
    axes[0].imshow(matrix, cmap='hot', aspect='auto')
    axes[0].set_xlabel('genomic region')
    axes[0].set_ylabel('genomic region')
    axes[0].set_title('Contact Matrix Heatmap')

    # Plot scatter plot
    n = len(labels)
    x = np.arange(n)
    y = np.zeros(n)
    axes[1].scatter(x, y, c=labels, cmap='Set1')
    axes[1].set_xlabel('genomic region')
    axes[1].set_ylabel('Cluster Label')
    axes[1].set_title('Cluster Labels Plot')

    # Adjust spacing between subplots
    fig.tight_layout()

    # Save combined image
    plt.savefig(filename)
    plt.close()
'''


# Parameters
n = 100  # Length of the vector
k = 10   # Number of entries per cluster (0's and 1's)
heatmap_filename = '/projects/b1198/epifluidlab/david/GSE63525/GM12878/Images/CheckerboardExample.png'  # Output filename
scatterplot_filename = '/projects/b1198/epifluidlab/david/GSE63525/GM12878/Images/Scatterplot.png'
combined_filename = '/projects/b1198/epifluidlab/david/GSE63525/GM12878/Images/combined_image.png'

# Create the cluster labels vector
labels = np.concatenate([np.zeros(k), np.ones(k)] * (n // (2 * k)))
labels2 = np.tile(np.repeat([0, 1], [15, 5]), n // 20)
labels3 = np.concatenate([np.repeat([0, 1], [k*i, k*(i+1)]) for i in range(n // k)])
labels4 = np.concatenate([np.zeros(k), np.ones(k), np.full(k, 2)] * (n // (3 * k)))


# Create the contact matrix
contact_matrix = create_contact_matrix(labels)

# Save the heatmap image
save_heatmap_image(contact_matrix, heatmap_filename)

# Save the scatterplot image
save_scatterplot_image(labels4, scatterplot_filename)

# Save the combined plots image
save_combined_image(contact_matrix, labels, combined_filename)






