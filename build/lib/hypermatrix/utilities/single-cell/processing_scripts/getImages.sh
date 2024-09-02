#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=1:48:00
#SBATCH --ntasks=5
#SBATCH --mem=90000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetImages.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JgetImages.out

data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"

python << EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import textwrap
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
from matplotlib.colors import Normalize
from matplotlib import cm

path = '$data_path'
resolutions = [1000000, 500000]

for resolution in resolutions:
    loaded_data = np.load(path + f'Images/dataForGraphs_res{resolution}.npy', allow_pickle=True).item()
    data_labels = loaded_data['data_labels']
    bin_averages_data_by_ch = loaded_data['bin_averages_data_by_ch']
    removed_rowIndices_by_ch = loaded_data['removed_rowIndices_by_ch']
    components_by_channel = loaded_data['components_by_channel']
    
    num_channels = 22
    channel_numbers = list(range(1, num_channels+1))
    sol_factors_downsampled_by_ch = {}
    components_by_channel = {ch: [] for ch in range(1, num_channels + 1)}
    data_ids = {'H3K9ac' : ['active', 'ENCFF128UVW'], 'H3K4me3': ['active','ENCFF154XCY'], 'H3K27me3':['repressive', 'ENCFF552FEL'], 'H3K9me3':['repressive', 'ENCFF533NIQ'], 'H4K20me1':['repressive','ENCFF646GIV'], 'DNAseI':['accessibility','ENCFF901GZH']}

    for ch in channel_numbers:
        mat_file = path + f'Workspaces/ch{ch}_res{resolution}_structedData_10Downsampled_pearsons_first3_500iterations.mat'
        data = scipy.io.loadmat(mat_file)
        sol_factors_downsampled_by_ch[ch] = data['sol_factors_downsampled']
        for component in range(1, 3):
            component_data = data['sol_factors_downsampled'][0, 1][0, 0][:, component - 1]
            components_by_channel[ch].append(component_data)

    def density_scatter(x, y, ax=None, sort=True, bins=20, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
        z = interpn((0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])), data, np.vstack([x, y]).T,
                method="splinef2d", bounds_error=False)
        z[np.where(np.isnan(z))] = 0.0
        if sort:
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
        sc = ax.scatter(x, y, c=z, **kwargs)
        norm = Normalize(vmin=np.min(z), vmax=np.max(z))
        cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax)
        cbar.ax.set_ylabel('Density', fontsize = 55)
        cbar.ax.tick_params(labelsize=55)
        return ax

    def create_scatter_plot(ax, filtered_y, filtered_binAveragesData, title):
        coefficients_deg1 = np.polyfit(filtered_y, filtered_binAveragesData, 1)
        best_fit_deg1 = np.poly1d(coefficients_deg1)
        x_values = np.linspace(min(filtered_y), max(filtered_y), 100)
        y_values_deg1 = best_fit_deg1(x_values)
        density_scatter(filtered_y, filtered_binAveragesData, ax=ax, sort=True, bins=20, cmap='viridis', s=300, alpha=0.7)
        ax.plot(x_values, y_values_deg1, color='magenta', linewidth=20, label='Best-Fit Line')
        pearson_corr, pearson_pval = pearsonr(filtered_y, filtered_binAveragesData)
        spearman_corr, spearman_pval = spearmanr(filtered_y, filtered_binAveragesData)
        ax.set_title(f'{title}\nPearson: {pearson_corr:.03f}, {pearson_pval:.03f}\nSpearman: {spearman_corr:.03f}, {spearman_pval:.03f}', fontsize=55)
        ax.set_xlabel('Euchromatin Factor', fontsize=55)
        ax.set_ylabel('Epigenetic Mark', fontsize=55)
        ax.tick_params(axis='both', labelsize=55)
        ax.legend(fontsize=45)

    def create_scatter_plot_aiden(ax, filtered_y, filtered_binAveragesData, title):
        coefficients_deg1 = np.polyfit(filtered_y, filtered_binAveragesData, 1)
        best_fit_deg1 = np.poly1d(coefficients_deg1)
        x_values = np.linspace(min(filtered_y), max(filtered_y), 100)
        y_values_deg1 = best_fit_deg1(x_values)
        density_scatter(filtered_y, filtered_binAveragesData, ax=ax, sort=True, bins=20, cmap='viridis', s=300, alpha=0.7)
        ax.plot(x_values, y_values_deg1, color='magenta', linewidth=20, label='Best-Fit Line')
        pearson_corr, pearson_pval = pearsonr(filtered_y, filtered_binAveragesData)
        spearman_corr, spearman_pval = spearmanr(filtered_y, filtered_binAveragesData)
        ax.set_title(f'{title}\nPearson: {pearson_corr:.03f}, {pearson_pval:.03f}\nSpearman: {spearman_corr:.03f}, {spearman_pval:.03f}', fontsize=55)
        ax.set_xlabel('Aiden Eigenvector', fontsize=55)
        ax.set_ylabel('Epigenetic Mark', fontsize=55)
        ax.tick_params(axis='both', labelsize=55)
        ax.legend(fontsize=45)

    #filter and concatenate
    filtered_data = {label: [] for label in data_labels}
    filtered_aiden = []
    filtered_y_all = []
    filtered_y_by_cat = {f'{label}_cat{i}': [] for label in data_labels for i in range(1, 5)}


    eigenvector_files = {
        ch: path + f'eigenvector/res{elementsPerBin}_ch{ch}_observed_NONE_NONE_eigenvector.txt'
        for ch in range(1, 23)
    }

    aiden_eigenvectors = {
        f'ch{ch}': load_and_replace_nan(file_path) for ch, file_path in eigenvector_files.items(
    )
    }



    for ch_num in channel_numbers:
        N_ch = sol_factors_downsampled_by_ch[ch_num][0, 1][0, 0][:, 0].size
        filtered_indices_ch = [i for i in range(N_ch-2) if i not in removed_rowIndices_by_ch[ch_num]]
        vector1 = components_by_channel[ch_num][0][filtered_indices_ch]
        vector2 = components_by_channel[ch_num][1][filtered_indices_ch]
        vector3 = bin_averages_data_by_ch[ch_num]['H3K9ac'][filtered_indices_ch]
        correlation1 = np.corrcoef(vector1, vector3)[0, 1]
        correlation2 = np.corrcoef(vector2, vector3)[0, 1]
        if correlation1 > correlation2:
            filtered_y = vector1
        else:
            filtered_y = vector2
        hypermatrix_mask = filtered_y > np.mean(filtered_y)
        filtered_y_all.extend(filtered_y)

        vector1 = aiden_eigenvectors[f'ch{ch_num}'][filtered_indices_ch]
        vector2 = -vector1
        correlation1 = np.corrcoef(vector1, vector3)[0, 1]
        correlation2 = np.corrcoef(vector2, vector3)[0, 1]
        if correlation1 > correlation2:
            filtered_aiden_eigenvectors = vector1
        else:
            filtered_aiden_eigenvectors = vector2
        aiden_mask = filtered_aiden_eigenvectors > 0
        filtered_aiden.extend(filtered_aiden_eigenvectors)
        for label in data_labels:
            filtered_data[label].extend(bin_averages_data_by_ch[ch_num][label][filtered_indices_ch])
        bin_values = np.array(filtered_indices_ch)
        cat1_values = bin_values[hypermatrix_mask & aiden_mask]
        cat2_values = bin_values[hypermatrix_mask & ~aiden_mask]
        cat3_values = bin_values[~hypermatrix_mask & aiden_mask]
        cat4_values = bin_values[~hypermatrix_mask & ~aiden_mask]

        for label in data_labels:
            filtered_y_by_cat[f'{label}_cat1'].extend(bin_averages_data_by_ch[ch_num][label][cat1_values])
            filtered_y_by_cat[f'{label}_cat2'].extend(bin_averages_data_by_ch[ch_num][label][cat2_values])
            filtered_y_by_cat[f'{label}_cat3'].extend(bin_averages_data_by_ch[ch_num][label][cat3_values])
            filtered_y_by_cat[f'{label}_cat4'].extend(bin_averages_data_by_ch[ch_num][label][cat4_values])

    # Hypermatrix Factor Figure    
    fig = plt.figure(figsize=(100, 55))
    for i in range(1,len(data_labels)+1):
        ax = fig.add_subplot(2, 3, i)
        title = f'{data_labels[i-1]} - {data_ids[data_labels[i-1]][0]}\n{data_ids[data_labels[i-1]][1]}'
        create_scatter_plot(ax, np.array(filtered_y_all), np.array(filtered_data[data_labels[i-1]]), title)

    fig.subplots_adjust(bottom=.25)
    fig.subplots_adjust(top=0.85)

    caption = f'These scatter plots show the purported euchromatin hypermatrix factor vs. epigenetic marks genome-wide. Each point (x, y) represents a genomic bin on an autosomal chromosome, where x is the value of the hypermatrix euchromatin factor at that bin and y is the value of the epigenetic factor at that bin. The hypermatrix factor was derived from 22 Hi-C contact matrices, one for each autosomal chromosome, at a resolution of {resolution} from data set GSE63525. The dataset is from GM12878 cells with only MAPQ > 30 reads and is labeled GSE63525_GM12878_insitu_primary+replicate_combined_30.hic. From this dataset, for each chromosome, the Pearson correlation matrix of the observed/expected matrix was computed using Straw and Juicer Tools. The matrix was translated by 1 to ensure non-negativity. Next, an n by n by 10 hypermatris was formed for each chromosome by downsampling, where n is the number of bins in that chromosome. Downsampling as done at 10 percent, 20 percent, and so on. Tensorlab\'s Structure Data Fusion sdf_nls algorithm was used to find an approximation of the optimal non-negative rank two decomposition with first two components symmetric. This produced two n by 1 vectors, referred to as the hypermatrix factors. The factor that was more correlated with H3K9ac was labeled the euchromatin factor. The euchromatin factor for each chromosome are on the x-axis of the plots. These plots demonstrate that the euchromatin factor indeed correlates with euchromatin at resolution {resolution}. Rows of the contact matrix with more than 90 percent zero entries were filtered out of the scatter plot. This was almost always the centromeric region.'
    caption_wrapped = textwrap.fill(caption, 200)
    fig.text(0.1, 0.03, caption_wrapped, ha='left', fontsize=55)

    fig.suptitle(f'Purported Euchromatin Factor and Epigenetic Marks Scatter Plots\nRes {resolution}, Chr1-22', fontsize=100, fontweight='bold')
    fig.subplots_adjust(hspace=1)
    plt.subplots_adjust(wspace=0.3)
    plt.tight_layout(rect=[0, 0.1, 1, 0.85])
    output_file_plot = path + f'Images/EpigeneticComparisonPlot_res{resolution}.png'
    plt.savefig(output_file_plot)
    plt.close()

    # Aidenlab Eigenvector Figure
    fig = plt.figure(figsize=(100, 55))
    for i in range(1,len(data_labels)+1):
        ax = fig.add_subplot(2, 3, i)
        title = f'{data_labels[i-1]} - {data_ids[data_labels[i-1]][0]}\n{data_ids[data_labels[i-1]][1]}'
        create_scatter_plot_aiden(ax, np.array(filtered_aiden), np.array(filtered_data[data_labels[i-1]]), title)
    fig.subplots_adjust(bottom=.25)
    fig.subplots_adjust(top=0.85)
    caption = f'These scatter plots show the correlation of the Aiden lab\'s Eigenvector command with epigenetic marks genome-wide. Each point (x, y) represents a genomic bin on an autosomal chromosome, where x is the value of the Eigenvector command at that bin and y is the value of the epigenetic factor at that bin. The inputs of the Eigenvector command were 22 autosomal chromosomes at a resolution of {resolution} from data set GSE63525. The dataset is from GM12878 cells with only MAPQ > 30 reads and is labeled GSE63525_GM12878_insitu_primary+replicate_combined_30.hic. The eigenvector outputed by the Aiden Eigenvector command was multiplied by negative one if that improved the vectors correlation with H3K9ac. These plots demonstrate that the hypermatrix euchromatin factor and the Aiden lab Eigenvector both correlate to epigenetic marks approximately the same amount. Rows of the contact matrix with more than 90 percent zero entries were filtered out of the scatter plot. This was almost always the centromeric region.'
    caption_wrapped = textwrap.fill(caption, 200)
    fig.text(0.1, 0.03, caption_wrapped, ha='left', fontsize=55)

    fig.suptitle(f'Aiden Eigenvector and Epigenetic Marks Scatter Plots\nRes {resolution}, Chr1-22', fontsize=100, fontweight='bold')
    fig.subplots_adjust(hspace=1)
    plt.subplots_adjust(wspace=0.3)
    plt.tight_layout(rect=[0, 0.1, 1, 0.85])
    output_file_plot = path + f'Images/EpigeneticAidenComparisonPlot_res{resolution}.png'
    plt.savefig(output_file_plot)
    plt.close()

    def plot_density_histograms(data_labels, filtered_data_by_cat, output_file_plot):
        fig, axs = plt.subplots(len(data_labels), 4, figsize=(15, 10))
        fig.suptitle('Density Histograms for Different Epigenetic Marks')

        # Calculate x_min and x_max for all categories and labels
        x_min = 0
        # Iterate over data_labels and plot each label in a separate subplot
        for i, label in enumerate(data_labels):
            x_max = max(max(filtered_y_by_cat[f'{label}_cat1']), max(filtered_y_by_cat[f'{label}_cat2']), max(filtered_y_by_cat[f'{label}_cat3']), max(filtered_y_by_cat[f'{label}_cat4']))
            # Category 1: Hypermatrix > Threshold, Aiden > 0
            axs[i, 0].hist(filtered_y_by_cat[f'{label}_cat1'], bins=10, density=False, alpha=0.7, label='Hypermatrix > Threshold, Aiden > 0', edgecolor='black')
            axs[i, 0].set_xlabel(f'{label} value')
            axs[i, 0].set_ylabel('Number of bins')
            axs[i, 0].legend(fontsize=6)
            axs[i, 0].set_xlim(x_min, x_max)

            # Category 2: Hypermatrix > Threshold, Aiden <= 0
            axs[i, 1].hist(filtered_y_by_cat[f'{label}_cat2'], bins=10, density=False, alpha=0.7, label='Hypermatrix > Threshold, Aiden <= 0', edgecolor='black')
            axs[i, 1].set_xlabel(f'{label} value')
            axs[i, 1].set_ylabel('Number of bins')
            axs[i, 1].legend(fontsize=6)
            axs[i, 1].set_xlim(x_min, x_max)

            # Category 3: Hypermatrix <= Threshold, Aiden > 0
            axs[i, 2].hist(filtered_y_by_cat[f'{label}_cat3'], bins=10, density=False, alpha=0.7, label='Hypermatrix <= Threshold, Aiden > 0', edgecolor='black')
            axs[i, 2].set_xlabel(f'{label} value')
            axs[i, 2].set_ylabel('Number of bins')
            axs[i, 2].legend(fontsize=6)
            axs[i, 2].set_xlim(x_min, x_max)

            # Category 4: Hypermatrix <= Threshold, Aiden <= 0
            axs[i, 3].hist(filtered_y_by_cat[f'{label}_cat4'], bins=10, density=False, alpha=0.7, label='Hypermatrix <= Threshold, Aiden <= 0', edgecolor='black')
            axs[i, 3].set_xlabel(f'{label} value')
            axs[i, 3].set_ylabel('Number of bins')
            axs[i, 3].legend(fontsize=6)
            axs[i, 3].set_xlim(x_min, x_max)

        caption = 'These density histograms show that the Aiden lab\'s Eigenvector and the Hypermatrix euchromatin factor are approximately of the same accuracy at predicting euchromatin epigenetic marks genome-wide. Each bin is a genomic bin on the 22 autosomal chromosomes at a resolution of {resolution} from data set GSE63525. The dataset is from GM12878 cells with only MAPQ > 30 reads and is labeled GSE63525_GM12878_insitu_primary+replicate_combined_30.hic. Rows of the contact matrix with more than 90 percent zero entries were filtered out of the scatter plot. This was almost always the centromeric region.'
        caption_wrapped = textwrap.fill(caption, 200)
        fig.text(0.01, 0.01, caption_wrapped, ha='left', fontsize=10)
        plt.tight_layout(rect=[0, 0.1, 1, 0.95])
        plt.savefig(output_file_plot)
        plt.close()

    # Create one figure with four plots for each label
    output_file_plot = path+'Images/HistogramPlot_res{resolution}.png'
    plot_density_histograms(data_labels, filtered_y_by_cat, output_file_plot)
EOF



python << Density

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import textwrap
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy.stats import ttest_ind, ks_2samp, kstwobign

resolutions = [1000000, 500000]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
chromosomes_as_integers = [int(ch) for ch in chromosomes]
data_labels = ['H3K9ac', 'H3K4me3', 'H3K27me3', 'H3K9me3', 'H4K20me1', 'DNAseI']
data_ids = {'H3K9ac' : ['active', 'ENCFF128UVW'], 'H3K4me3': ['active','ENCFF154XCY'], 'H3K27me3':['repressive', 'ENCFF552FEL'], 'H3K9me3':['repressive', 'ENCFF533NIQ'], 'H4K20me1':['repressive','ENCFF646GIV'], 'DNAseI':['accessibility','ENCFF901GZH']}
path = '/data/epifluidlab/david/GSE63525version2/GM12878/'
num_channels = len(chromosomes)
channel_numbers = list(range(1, num_channels+1))

# Function to create density plot
def create_density_plot_3clusters(allData, label, id, resolution):
    fig, axs = plt.subplots(2, 3, figsize=(15, 15))
    fig.suptitle(f'Density Histograms for {label}\n{id[0]}, {id[1]}, res {resolution:,}', fontsize = 35)
    # Compute means and standard deviations for each distribution
    mean_cat1 = np.mean(allData[f'{label}_cat1'])
    std_cat1 = np.std(allData[f'{label}_cat1'])
    mean_cat2 = np.mean(allData[f'{label}_cat2'])
    std_cat2 = np.std(allData[f'{label}_cat2'])
    mean_cat3 = np.mean(allData[f'{label}_cat3'])
    std_cat3 = np.std(allData[f'{label}_cat3'])
    mean_cat4 = np.mean(allData[f'{label}_cat4'])
    std_cat4 = np.std(allData[f'{label}_cat4'])
    mean_cat5 = np.mean(allData[f'{label}_cat5'])
    std_cat5 = np.std(allData[f'{label}_cat5'])
    mean_cat6 = np.mean(allData[f'{label}_cat6'])
    std_cat6 = np.std(allData[f'{label}_cat6'])

    # Calculate total number of bins for each category
    total_bins_cat1 = len(allData[f'{label}_cat1'])
    total_bins_cat2 = len(allData[f'{label}_cat2'])
    total_bins_cat3 = len(allData[f'{label}_cat3'])
    total_bins_cat4 = len(allData[f'{label}_cat4'])
    total_bins_cat5 = len(allData[f'{label}_cat5'])
    total_bins_cat6 = len(allData[f'{label}_cat6'])


    # Significance level (alpha)
    alpha = 0.05

    # Number of simulation trials
    num_trials = 10000  # Adjust the number of trials as needed

    # Sample sizes for the two datasets
    n1 = len(allData[f'{label}_cat5'])  # Use the appropriate sample size
    n2 = len(allData[f'{label}_cat6'])  # Use the appropriate sample size

    # Initialize an array to store K-S statistics from simulations
    ks_statistics = np.zeros(num_trials)

    # Perform simulations
    for i in range(num_trials):
        # Generate random samples from suitable distributions
        sample1 = np.random.normal(mean_cat5, std_cat5, n1)  # Use appropriate distributions and parameters
        sample2 = np.random.normal(mean_cat6, std_cat6, n2)  # Use appropriate distributions and parameters

        # Calculate K-S statistic for the simulated samples
        ks_statistic = ks_2samp(sample1, sample2).statistic
        ks_statistics[i] = ks_statistic

    # Calculate the empirical critical value based on simulations
    critical_value_empirical = np.percentile(ks_statistics, (1 - alpha / 2) * 100)

    # Calculate the critical value
    critical_value = kstwobign.ppf(1 - alpha / 2) * np.sqrt((n1 + n2) / (n1 * n2))

    x_min = 0

    # Gather all the sequences into a list
    sequences = [
        allData[f'{label}_cat1'],
        allData[f'{label}_cat2'],
        allData[f'{label}_cat3'],
        allData[f'{label}_cat4'],
        allData[f'{label}_cat5'],
        allData[f'{label}_cat6']
    ]
    # Filter out empty sequences and get the maximum value from each non-empty sequence
    max_values = [max(seq) for seq in sequences if seq]
    # Check if there are any non-empty sequences
    if max_values:
        x_max = max(max_values)
    else:
        x_max = 10
        print("No data in the sequences.")

    axs[0, 0].hist(allData[f'{label}_cat1'], bins=10, density=False, alpha=0.7, label='Hypermatrix Euchromatin, Aiden Euchromatin', edgecolor='black')
    axs[0, 0].set_xlabel(f'{label} value')
    axs[0, 0].set_ylabel('Number of bins')
    axs[0, 0].legend(fontsize=10)
    axs[0, 0].set_xlim(x_min, x_max)
    axs[0, 0].set_title(f'Mean: {mean_cat1:.2f}, Std: {std_cat1:.2f}\nTotal Bins: {total_bins_cat1}', fontsize = 15)
    axs[0, 1].hist(allData[f'{label}_cat2'], bins=10, density=False, alpha=0.7, label='Hypermatrix Euchromatin, Aiden Heterochromatin', edgecolor='black')
    axs[0, 1].set_xlabel(f'{label} value')
    axs[0, 1].set_ylabel('Number of bins')
    axs[0, 1].legend(fontsize=10)
    axs[0, 1].set_xlim(x_min, x_max)
    axs[0, 1].set_title(f'Mean: {mean_cat2:.2f}, Std: {std_cat2:.2f}\nTotal Bins: {total_bins_cat2}', fontsize = 15)
    axs[1, 0].hist(allData[f'{label}_cat3'], bins=10, density=False, alpha=0.7, label='Hypermatrix Heterochromatin, Aiden Euchromatin', edgecolor='black')
    axs[1, 0].set_xlabel(f'{label} value')
    axs[1, 0].set_ylabel('Number of bins')
    axs[1, 0].legend(fontsize=10)
    axs[1, 0].set_xlim(x_min, x_max)
    axs[1, 0].set_title(f'Mean: {mean_cat3:.2f}, Std: {std_cat3:.2f}\nTotal Bins: {total_bins_cat3}', fontsize = 15)
    axs[1, 1].hist(allData[f'{label}_cat4'], bins=10, density=False, alpha=0.7, label='Hypermatrix Heterochromatin, Aiden Heterochromatin', edgecolor='black')
    axs[1, 1].set_xlabel(f'{label} value')
    axs[1, 1].set_ylabel('Number of bins')
    axs[1, 1].legend(fontsize=10)
    axs[1, 1].set_xlim(x_min, x_max)
    axs[1, 1].set_title(f'Mean: {mean_cat4:.2f}, Std: {std_cat4:.2f}\nTotal Bins: {total_bins_cat4}', fontsize = 15)
    axs[0, 2].hist(allData[f'{label}_cat5'], bins=10, density=False, alpha=0.7, label='Hypermatrix Neither, Aiden Euchromatin', edgecolor='black')
    axs[0, 2].set_xlabel(f'{label} value')
    axs[0, 2].set_ylabel('Number of bins')
    axs[0, 2].legend(fontsize=10)
    axs[0, 2].set_xlim(x_min, x_max)
    axs[0, 2].set_title(f'Mean: {mean_cat5:.2f}, Std: {std_cat5:.2f}\nTotal Bins: {total_bins_cat5}', fontsize = 15)
    axs[1, 2].hist(allData[f'{label}_cat6'], bins=10, density=False, alpha=0.7, label='Hypermatrix Neither, Aiden Heterochromatin', edgecolor='black')
    axs[1, 2].set_xlabel(f'{label} value')
    axs[1, 2].set_ylabel('Number of bins')
    axs[1, 2].legend(fontsize=10)
    axs[1, 2].set_xlim(x_min, x_max)
    axs[1, 2].set_title(f'Mean: {mean_cat6:.2f}, Std: {std_cat6:.2f}\nTotal Bins: {total_bins_cat6}', fontsize = 15)

    p_value_cat1_vs_cat4 = ttest_ind(allData[f'{label}_cat1'], allData[f'{label}_cat4'], equal_var=False).pvalue

    # Run K-S test for categories 5 and 6
    ks_statistic, ks_p_value = ks_2samp(allData[f'{label}_cat5'], allData[f'{label}_cat6'])
    caption = f'These density histograms show the comparison of a 3 Clustering with two hypermatrix factors with Aiden Eigenvector 2 clustering for chr 1-22 with respect to {label} with a resolution of {resolution}. Centromeric bins are not included in the clustering.'
    caption += f' The p-value for the difference in mean between Histogram in (1,1) and Histogram in (2,2) is: {p_value_cat1_vs_cat4:.4f}.'
    caption += f' The K-S statistic for Histogram in (1,3) and Histogram in (2,3) is: {ks_statistic:.4f} with p-value: {ks_p_value:.4f}. A large p-value is expected due to small sample size.'
    caption += f' For the K-S statistic, the critical value: {critical_value:.3f}, and empirical critical value: {critical_value_empirical:.3f}.'
    caption += f' Each bin is a genomic bin on the 22 autosomal chromosomes at a resolution of {resolution} from data set GSE63525. The dataset is from GM12878 cells with only MAPQ > 30 reads and is labeled GSE63525_GM12878_insitu_primary+replicate_combined_30.hic. Rows of the contact matrix with more than 90 percent zero entries were filtered out of the scatter plot. This was almost always the centromeric region.'
    caption_wrapped = textwrap.fill(caption, 98)
    fig.text(0.01, 0.01, caption_wrapped, ha='left', fontsize=20)
    plt.tight_layout(rect=[0, 0.25, 1, 0.97])
    output_file_plot = path+f'Images/HistogramPlot_res{resolution}_allChromosomes_{label}_3clustering.png'
    plt.savefig(output_file_plot)
    plt.close()
for resolution in resolutions:
    allData = {f'{label}_cat{i}': [] for label in data_labels for i in range(1, 7)}
    loaded_data = np.load(path + f'Images/dataForGraphs_res{resolution}.npy', allow_pickle=True).item()
    bin_averages_data_by_ch = loaded_data['bin_averages_data_by_ch']
    removed_rowIndices_by_ch = loaded_data['removed_rowIndices_by_ch']
    components_by_channel = loaded_data['components_by_channel']
    aiden_eigenvectors = loaded_data['aiden_eigenvectors']
    sol_factors_downsampled_by_ch = {}
    for ch in chromosomes_as_integers:
        N_ch = components_by_channel[ch][0].size
        filtered_indices_ch = [i for i in range(N_ch-1) if i not in removed_rowIndices_by_ch[ch]]
        vector1 = components_by_channel[ch][0][filtered_indices_ch]
        vector2 = components_by_channel[ch][1][filtered_indices_ch]
        vector3 = bin_averages_data_by_ch[ch]['H3K9ac'][filtered_indices_ch]
        correlation1 = np.corrcoef(vector1, vector3)[0, 1]
        correlation2 = np.corrcoef(vector2, vector3)[0, 1]
        if correlation1 > correlation2:
            euchromatin_factor = components_by_channel[ch][0]
            heterochromatin_factor = components_by_channel[ch][1]
        else:
            euchromatin_factor = components_by_channel[ch][1]
            heterochromatin_factor = components_by_channel[ch][0]
        euchromatin_mask = np.logical_and(euchromatin_factor > np.mean(euchromatin_factor[filtered_indices_ch]), heterochromatin_factor < np.mean(heterochromatin_factor[filtered_indices_ch]))
        heterochromatin_mask = np.logical_and(euchromatin_factor < np.mean(euchromatin_factor[filtered_indices_ch]), heterochromatin_factor > np.mean(heterochromatin_factor[filtered_indices_ch]))
        neither_mask = ~(euchromatin_mask | heterochromatin_mask)
        vector1 = aiden_eigenvectors[f'ch{ch}'][filtered_indices_ch]
        vector2 = -vector1
        correlation1 = np.corrcoef(vector1, vector3)[0, 1]
        correlation2 = np.corrcoef(vector2, vector3)[0, 1]
        if correlation1 > correlation2:
            aiden_factor = aiden_eigenvectors[f'ch{ch}']
        else:
            aiden_factor = -aiden_eigenvectors[f'ch{ch}']
        aiden_euchromatin = aiden_factor > 0
        aiden_heterochromatin = aiden_factor < 0
        bins = np.array(filtered_indices_ch)
        cat1_bin_indices = bins[euchromatin_mask[filtered_indices_ch] & aiden_euchromatin[filtered_indices_ch]]
        cat2_bin_indices = bins[euchromatin_mask[filtered_indices_ch] & ~aiden_euchromatin[filtered_indices_ch]]
        cat3_bin_indices = bins[heterochromatin_mask[filtered_indices_ch] & aiden_euchromatin[filtered_indices_ch]]
        cat4_bin_indices = bins[heterochromatin_mask[filtered_indices_ch] & ~aiden_euchromatin[filtered_indices_ch]]
        cat5_bin_indices = bins[neither_mask[filtered_indices_ch] & aiden_euchromatin[filtered_indices_ch]]
        cat6_bin_indices = bins[neither_mask[filtered_indices_ch] & ~aiden_euchromatin[filtered_indices_ch]]
        for label in data_labels:
            allData[f'{label}_cat1'].extend(bin_averages_data_by_ch[ch][label][cat1_bin_indices])
            allData[f'{label}_cat2'].extend(bin_averages_data_by_ch[ch][label][cat2_bin_indices])
            allData[f'{label}_cat3'].extend(bin_averages_data_by_ch[ch][label][cat3_bin_indices])
            allData[f'{label}_cat4'].extend(bin_averages_data_by_ch[ch][label][cat4_bin_indices])
            allData[f'{label}_cat5'].extend(bin_averages_data_by_ch[ch][label][cat5_bin_indices])
            allData[f'{label}_cat6'].extend(bin_averages_data_by_ch[ch][label][cat6_bin_indices])
    for label in data_labels:
        id = data_ids[label]
        create_density_plot_3clusters(allData, label, id, resolution)

Density










