#!/bin/bash
#for interactive session: srun -N 1 -n 1 --account=b1042 --mem=50000 --partition=genomics --time=1:00:00 --pty bash -l


#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=1:10:00
#SBATCH --ntasks=5
#SBATCH --mem=50000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%J.insulationImages.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%J.insulationImages.out

data_path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"

resolutions=(1000000)
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")
#chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
data_types=("observed")
genomeID="hg19"

resolutions_list=$(printf "'%s', " "${resolutions[@]}")
chromosomes_list=$(printf "'%s', " "${chromosomes[@]}")
data_types_list=$(printf "'%s', " "${data_types[@]}")

resolutions_list=${resolutions_list%, }
chromosomes_list=${chromosomes_list%, }
data_types_list=${data_types_list%, }



python << getData
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import textwrap
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
from matplotlib.colors import Normalize
from matplotlib import cm
import pyBigWig
import h5py
import math

path = '$data_path'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
insulation_windows = ["375kb", "250kb", "175kb", "125kb", "75kb"] #for res25,000
#insulation_windows = ["7.5mb", "5mb", "3.5mb", "2.5mb", "1.5mb"] #for res500,000
sol_factors_downsampled_by_ch = {}
removed_rowIndices_by_ch = {}
insulation_by_window_ch = {window : {} for window in insulation_windows}

def calculate_bin_averages(data, elements_per_bin):
    num_bins = math.ceil(len(data)/elements_per_bin)
    bin_averages = np.zeros(num_bins)
    for i in range(num_bins):
        start_index = i * elements_per_bin
        end_index = min((i + 1) * elements_per_bin, len(data))
        bin_data = data[start_index:end_index]
        if len(bin_data) > 0:
            bin_averages[i] = np.mean(bin_data)
        else:
            bin_averages[i] = 0
    return np.nan_to_num(bin_averages, nan=0.0)

#store the matrices, solution factors, and removed rows in dictionaries
for resolution in resolutions:
    resolution = int(resolution)   
    for ch in chromosomes:
        mat_file = path + f'Workspaces/individual/ch{ch}_res{resolution}_structedData_10Downsampled_pearsons_first3_500iterations.mat'
        data = scipy.io.loadmat(mat_file)
        sol_factors_downsampled_by_ch[ch] = data['sol_factors_downsampled']
        for window in insulation_windows:
            bigwig_file = path + f'insulation/insulation_{resolution}_ch{ch}_observed_KR_{window}.bigwig'
            bw = pyBigWig.open(bigwig_file)
            start = 1
            end = bw.chroms()[ch]
            insulation_by_window_ch[window][ch] = calculate_bin_averages(bw.values(ch, start, end), resolution)

    # Gather all variables you want to save into a single dictionary
    data_to_save = {
        'insulation_by_window_ch': insulation_by_window_ch,
        'sol_factors_downsampled_by_ch': sol_factors_downsampled_by_ch
        }
    # Save the data to an npy file
    np.save(path + f'Images/dataForInsulation_res{resolution}.npy', data_to_save)
getData


<<comment
python<<EOF
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import textwrap
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
from matplotlib.colors import Normalize
from matplotlib import cm
import pyBigWig
import h5py
import math
from scipy.stats import t
from scipy.signal import convolve

path = '$data_path'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
insulation_windows = ["375kb", "250kb", "175kb", "125kb", "75kb"]

# Define a list of smoothing factors to try for hypermatrix vector and insulation vector
hypermatrix_smoothing_factors = [1, 10, 100, 510, 600, 700, 1000, 2000]  # Adjust as needed
insulation_vector_smoothing_factors = [1, 10, 100, 510, 500, 600, 700, 1000, 2000]  # Adjust as needed

def smooth_vector(vector, smoothing_factor):
    kernel = np.ones(smoothing_factor) / smoothing_factor
    smoothed_vector = convolve(vector, kernel, mode='same')
    return smoothed_vector

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
        cbar.ax.set_ylabel('Density', fontsize = 10)
        cbar.ax.tick_params(labelsize=10)
        return ax

def create_scatter_plot(ax, filtered_y, filtered_binAveragesData, title):
        coefficients_deg1 = np.polyfit(filtered_y, filtered_binAveragesData, 1)
        best_fit_deg1 = np.poly1d(coefficients_deg1)
        x_values = np.linspace(min(filtered_y), max(filtered_y), 100)
        y_values_deg1 = best_fit_deg1(x_values)
        density_scatter(filtered_y, filtered_binAveragesData, ax=ax, sort=True, bins=20, cmap='viridis', s=50, alpha=0.7)
        ax.plot(x_values, y_values_deg1, color='magenta', linewidth=5, label='Best-Fit Line')
        pearson_corr, pearson_pval = pearsonr(filtered_y, filtered_binAveragesData)
        spearman_corr, spearman_pval = spearmanr(filtered_y, filtered_binAveragesData)
        ax.set_title(f'{title}\nPearson: {pearson_corr:.03f}, {pearson_pval:.03f}\nSpearman: {spearman_corr:.03f}, {spearman_pval:.03f}', fontsize=10)
        ax.set_xlabel('Hypermatrix Factor', fontsize=15)
        ax.set_ylabel('Insulation Score', fontsize=15)
        ax.tick_params(axis='both', labelsize=15)
        ax.legend(fontsize=10)

for resolution in resolutions:
    loaded_data = np.load(path + f'Images/dataForInsulation_res{resolution}_oe_KR_log2.npy', allow_pickle=True).item()
    insulation_by_window_ch = loaded_data['insulation_by_window_ch']
    sol_factors_downsampled_by_ch = loaded_data['sol_factors_downsampled_by_ch']
    matrices = {}
    removed_rowIndices_by_ch = {}

    for ch in chromosomes:
        vector1 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 0]
        vector2 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 1]
        vector3 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 2]
        vector4 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 3]
        #vector5 = sol_factors_downsampled_by_ch[ch][0, 4][0, 0][:, 4]


        #normalization_factor_two = np.sqrt(sol_factors_downsampled_by_ch[ch][0, 1][0, 0][0, 1] / sol_factors_downsampled_by_ch[ch][0, 1][0, 0][0, 0])
        #vector2 = normalization_factor_two * vector2

        #vector6 = np.exp((vector1 - vector2)**2)

        # Initialize variables to store the most correlated vectors and their correlation values
        max_corr_vector_pair = None
        max_corr_value = -1
        max_corr_window = None
        max_corr_hypermatrix = None
        max_hypermatrix_smoothing_factor = None
        max_insulation_smoothing_factor = None

        removed_rows_path = path + f'Workspaces/individual/ch{ch}_res{resolution}_oe_KR_removedRows.mat'    
        data = scipy.io.loadmat(removed_rows_path)
        removed_indices = data['removed_rows_indices']
        # Convert the indices to a NumPy array
        result_array = np.array(removed_indices)
        removed_rowIndices_by_ch[ch] = result_array
        N_ch = sol_factors_downsampled_by_ch[ch][0, 1][0, 0][:, 0].size
        filtered_indices_ch = [i for i in range(N_ch-1) if i not in removed_rowIndices_by_ch[ch]]

        #vector_mean = np.mean(vector6)

        # Take the base-2 logarithm of the result_vector
        #log2_result_vector = np.log2(vector6)

        # Loop through hypermatrix smoothing factors
        for hypermatrix_smoothing_factor in hypermatrix_smoothing_factors:
            # Loop through insulation vector smoothing factors
            for insulation_vector_smoothing_factor in insulation_vector_smoothing_factors:
                for idx, vector in enumerate([vector1, vector2, vector3, vector4]):
                    for window in insulation_windows:
                        # Extract the insulation vector from the second set
                        insulation_vector = insulation_by_window_ch[window][ch]

                        # Define the threshold below which values will be set
                        threshold = -5  # You can adjust the threshold as needed

                        # Set values below the threshold to be equal to the threshold
                        insulation_vector[insulation_vector < threshold] = threshold
                        

                        # Define the threshold below which values will be set
                        positive_threshold = 40  # You can adjust the threshold as needed

                        # Set values below the threshold to be equal to the threshold
                        vector[vector > positive_threshold] = positive_threshold

                        # Smooth both the hypermatrix vector and insulation vector
                        smoothed_hypermatrix_vector = smooth_vector(vector, hypermatrix_smoothing_factor)
                        smoothed_insulation_vector = smooth_vector(insulation_vector, insulation_vector_smoothing_factor)

                        # Calculate Pearson correlation coefficient and p-value for smoothed vectors
                        corr_value, p_value = pearsonr(smoothed_hypermatrix_vector[filtered_indices_ch], smoothed_insulation_vector[filtered_indices_ch])

                        # Update best correlation results if a better correlation is found
                        if corr_value > max_corr_value:
                            max_corr_value = corr_value
                            max_corr_window = window
                            max_corr_hypermatrix = idx + 1
                            max_corr_vector_pair = (vector, insulation_vector)
                            max_hypermatrix_smoothing_factor = hypermatrix_smoothing_factor
                            max_insulation_smoothing_factor = insulation_vector_smoothing_factor

        #loaded_data = np.load(path + f'Images/dataForGraphs_res{resolution}.npy', allow_pickle=True).item()
        #bin_averages_data_by_ch = loaded_data['bin_averages_data_by_ch']
        #factor1 = vector1[filtered_indices_ch]
        #factor2 = vector2[filtered_indices_ch]
        #vector3 = bin_averages_data_by_ch[ch_num]['H3K9ac'][filtered_indices_ch]
        #correlation1 = np.corrcoef(vector1, vector3)[0, 1]
        #correlation2 = np.corrcoef(vector2, vector3)[0, 1]
        #if correlation1 > correlation2:
        #    filtered_y = vector1
        #else:
        #    filtered_y = vector2

        mat_data = path + f'Workspaces/individual/ch{ch}_res{resolution}_oe_KR_matrix.mat'
        file = h5py.File(mat_data, 'r')
        matrix = file['matrix'][:]
        matrices[ch] = matrix

        # Plot the two most correlated vectors
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        # Scatter plot on the first panel
        #axs[0].scatter(most_corr_vectors[0], most_corr_vectors[1], label=f"Correlation: {max_corr_value:.2f}, p-value: {p_value:.2f}")
        #axs[0].set_title("Scatter Plot")
        axs[0].grid(True)
        title = "Scatter Plot"
        create_scatter_plot(axs[0], max_corr_vector_pair[0], max_corr_vector_pair[1], title)
        axs[0].set_xlabel("Hypermatrix Factor" + f"{max_corr_hypermatrix}")
        axs[0].set_ylabel("Insulation Score" + f" ({max_corr_window})")  

        # Additional plots on the second panel
        x_data = list(range(1, len(max_corr_vector_pair[1]) + 1))
        axs[1].plot(x_data[2800:3000], max_corr_vector_pair[0][2800:3000], label="Hypermatrix Factor" + f' {max_corr_hypermatrix}')
        axs[1].plot(x_data[2800:3000], max_corr_vector_pair[1][2800:3000], label="Insulation Score")
        axs[1].set_title("Factor Plots")
        axs[1].set_xlabel("Genomic Bin")
        axs[1].set_ylabel("Value")
        axs[1].legend()
        axs[1].grid(True)

        # heatmap for third  panel
        #vmin = 0
        #vmax = 5
        #heatmap = axs[2].imshow(matrix, cmap='viridis', vmin=vmin, vmax=vmax)
        #axs[2].set_title("Heatmap " + f"res{resolution} " + f"ch{ch}")
        #axs[2].set_xlabel("Genomic Bin")
        #axs[2].set_ylabel("Genomic Bin")
        #axs[2].grid(True)

        # Add a colorbar for reference
        #colorbar = plt.colorbar(heatmap, ax=axs, shrink=0.7)
        #colorbar.set_label("Color Scale Label")

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.3)

        caption_text = f"This figure shows the correlation with the centromeres removed of the hypermatrix factor\nfrom downsampled oe KR log2 data and FAN-C insulation score from observed KR data.\nHypermatrix Factor chosen from two most correlated and most correlated window size.\nHypermatrix Smoothing Factor: {max_hypermatrix_smoothing_factor}\nInsulating Vector Smoothing Factor: {max_insulation_smoothing_factor}"
        plt.figtext(0.5, 0.05, caption_text, ha="center", fontsize=12)
        plt.suptitle(f"Res {resolution} ch {ch}", fontsize=20)

        # Save the figure with both panels to a file
        plt.savefig(path + f'Images/res_{resolution}_ch{ch}_correlation_plot_log2.png')
        # Close the figure 
        plt.close(fig)
EOF




python<<wholeGenome
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import textwrap
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
from matplotlib.colors import Normalize
from matplotlib import cm
import pyBigWig
import h5py
import math
from scipy.stats import t
from scipy.signal import convolve

path = '$data_path'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
insulation_windows = ["375kb", "250kb", "175kb", "125kb", "75kb"]


def smooth_vector(vector, smoothing_factor):
    kernel = np.ones(smoothing_factor) / smoothing_factor
    smoothed_vector = convolve(vector, kernel, mode='same')
    return smoothed_vector

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
        cbar.ax.set_ylabel('Density', fontsize = 10)
        cbar.ax.tick_params(labelsize=10)
        return ax

def create_scatter_plot(ax, filtered_y, filtered_binAveragesData, title):
        coefficients_deg1 = np.polyfit(filtered_y, filtered_binAveragesData, 1)
        best_fit_deg1 = np.poly1d(coefficients_deg1)
        x_values = np.linspace(min(filtered_y), max(filtered_y), 100)
        y_values_deg1 = best_fit_deg1(x_values)
        density_scatter(filtered_y, filtered_binAveragesData, ax=ax, sort=True, bins=20, cmap='viridis', s=50, alpha=0.7)
        ax.plot(x_values, y_values_deg1, color='magenta', linewidth=5, label='Best-Fit Line')
        ax.set_title(f'{title}', fontsize=15)
        ax.set_xlabel('Hypermatrix Factor', fontsize=15)
        ax.set_ylabel('Insulation Score', fontsize=15)
        ax.tick_params(axis='both', labelsize=15)
        ax.legend(fontsize=10)

removed_rowIndices_by_ch = {}
matrices = {}

for resolution in resolutions:
    loaded_data = np.load(path + f'Images/dataForInsulation_res{resolution}_oe_KR_log2.npy', allow_pickle=True).item()
    insulation_by_window_ch = loaded_data['insulation_by_window_ch']
    sol_factors_downsampled_by_ch = loaded_data['sol_factors_downsampled_by_ch']
    
    for window in insulation_windows:
        filtered_hypermatrix_all = []
        filtered_insulation_all = []
        for ch in chromosomes:
                vector1 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 0]
                vector2 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 1]
                vector3 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 2]
                vector4 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 3]
                #vector5 = sol_factors_downsampled_by_ch[ch][0, 4][0, 0][:, 4]

                # Initialize variables to store the most correlated vectors and their correlation values
                most_corr_vectors = None
                max_corr_value = -1
                max_corr_window = None
                max_corr_hypermatrix = None

                removed_rows_path = path + f'Workspaces/individual/ch{ch}_res{resolution}_oe_KR_removedRows.mat'    
                data = scipy.io.loadmat(removed_rows_path)
                removed_indices = data['removed_rows_indices']
                #Convert the indices to a NumPy array
                result_array = np.array(removed_indices)
                removed_rowIndices_by_ch[ch] = result_array
                N_ch = sol_factors_downsampled_by_ch[ch][0, 1][0, 0][:, 0].size
                filtered_indices_ch = [i for i in range(N_ch-1) if i not in removed_rowIndices_by_ch[ch]]

                for idx, vector in enumerate([vector1, vector2, vector3, vector4]):
                    # Extract the insulation vector from the second set
                    insulation_vector = insulation_by_window_ch[window][ch]
                
                    # Define the threshold below which values will be set
                    threshold = -5  # You can adjust the threshold as needed

                    # Set values below the threshold to be equal to the threshold
                    insulation_vector[insulation_vector < threshold] = threshold

                    # Define the threshold below which values will be set
                    positive_threshold = 40  # You can adjust the threshold as needed

                    # Set values below the threshold to be equal to the threshold
                    vector[vector > positive_threshold] = threshold

                    # Smooth both the hypermatrix vector and insulation vector
                    smoothed_hypermatrix_vector = smooth_vector(vector, 100)
                    smoothed_insulation_vector = smooth_vector(insulation_vector, 100)

                    # Calculate Pearson correlation coefficient and p-value for smoothed vectors
                    corr_value, p_value = pearsonr(smoothed_hypermatrix_vector[filtered_indices_ch], smoothed_insulation_vector[filtered_indices_ch])

                    if corr_value > max_corr_value:
                        most_corr_vectors = (vector, insulation_vector)
                        max_corr_value = corr_value
                        max_corr_hypermatrix = idx+1

                filtered_hypermatrix_all.extend(most_corr_vectors[0][filtered_indices_ch])
                filtered_insulation_all.extend(most_corr_vectors[1][filtered_indices_ch])

        corr_value_all, p_value_all = pearsonr(filtered_hypermatrix_all, filtered_insulation_all)

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        title = f'Genome-wide scatter plot for window {window}'
        create_scatter_plot(ax, np.array(filtered_hypermatrix_all), np.array(filtered_insulation_all), title)
        plt.subplots_adjust(bottom=0.3)
        caption_text = f"This figure shows the genome-wide correlation with the centromeres removed\n of the hypermatrix factor from downsampled \noe KR log2 data and FAN-C insulation score from observed KR data.\nResolution {resolution} and window size {window} \nCorrelation: {corr_value_all:.2f}\nP-Value: {p_value_all:.2f}"
        plt.figtext(0.5, 0.1, caption_text, ha="center", fontsize=12)

        # Save the figure with both panels to a file
        plt.savefig(path + f'Images/res_{resolution}_window{window}_correlation_plot_genome_wide.png')
        # Close the figure 
        plt.close(fig)
wholeGenome





python<<correlatedSections
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import textwrap
from scipy.stats import pearsonr, spearmanr
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interpn
from matplotlib.colors import Normalize
from matplotlib import cm
import pyBigWig
import h5py
import math
from scipy.stats import t
from scipy.signal import correlate2d
import heapq
from matplotlib.gridspec import GridSpec

path = '$data_path'
resolutions = [$resolutions_list]
chromosomes = [$chromosomes_list]
insulation_windows = ["375kb", "250kb", "175kb", "125kb", "75kb"]
removed_rowIndices_by_ch = {}

def find_most_least_correlated_sections(vector1, vector2, window_size):
    n = len(vector1)
    if n != len(vector2):
        raise ValueError("Input vectors must have the same length")
    # Create a list to store the correlation values and their positions
    correlations = []
    # Iterate through the vectors with a sliding window
    for i in range(n - window_size + 1, window_size):
        section1 = vector1[i:i + window_size]
        section2 = vector2[i:i + window_size]
        # Calculate the correlation between the two sections
        correlation = np.corrcoef(section1, section2)[0, 1]
        # Store the correlation value and its position
        correlations.append((correlation, i))
    # Find the top 10 most correlated sections
    top_10_most_correlated = heapq.nlargest(5, correlations, key=lambda x: x[0])
    # Find the top 10 least correlated sections (negative correlation)
    top_10_least_correlated = heapq.nsmallest(5, correlations, key=lambda x: x[0])
    return top_10_most_correlated, top_10_least_correlated


for resolution in resolutions:
    resolution = int(resolution)
    insulation_window = "175kb"
    
    loaded_data = np.load(path + f'Images/dataForInsulation_res{resolution}_oe_KR_log2.npy', allow_pickle=True).item()
    insulation_by_window_ch = loaded_data['insulation_by_window_ch']
    sol_factors_downsampled_by_ch = loaded_data['sol_factors_downsampled_by_ch']

    for ch in chromosomes:
        vector1 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 0]
        vector2 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 1]
        vector3 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 2]
        vector4 = sol_factors_downsampled_by_ch[ch][0, 3][0, 0][:, 3]
        #vector5 = sol_factors_downsampled_by_ch[ch][0, 4][0, 0][:, 4]

        # Initialize variables to store the most correlated vectors and their correlation values
        most_corr_vectors = None
        max_corr_value = -1
        max_corr_window = None
        max_corr_hypermatrix = None

        removed_rows_path = path + f'Workspaces/individual/ch{ch}_res{resolution}_oe_KR_removedRows.mat'    
        data = scipy.io.loadmat(removed_rows_path)
        removed_indices = data['removed_rows_indices']
        #Convert the indices to a NumPy array
        result_array = np.array(removed_indices)
        removed_rowIndices_by_ch[ch] = result_array
        N_ch = sol_factors_downsampled_by_ch[ch][0, 1][0, 0][:, 0].size
        filtered_indices_ch = [i for i in range(N_ch-1) if i not in removed_rowIndices_by_ch[ch]]

        for idx, vector in enumerate([vector1, vector2, vector3, vector4]):
            # Extract the insulation vector from the second set
            insulation_vector = insulation_by_window_ch[insulation_window][ch]
                
            # Define the threshold below which values will be set
            threshold = -10  # You can adjust the threshold as needed

            # Set values below the threshold to be equal to the threshold
            insulation_vector[insulation_vector < threshold] = threshold

            # Calculate Pearson correlation coefficient and p-value
            corr_value, p_value = pearsonr(vector[filtered_indices_ch], insulation_vector[filtered_indices_ch])

            if corr_value > max_corr_value:
                most_corr_vectors = (vector, insulation_vector)
                max_corr_value = corr_value
                max_corr_hypermatrix = idx+1

        section_size = 1000;
        most_correlated, least_correlated = find_most_least_correlated_sections(most_corr_vectors[0], most_corr_vectors[1], section_size)

        # Create a 2x5 grid of subplots for the plots
        fig = plt.figure(figsize=(12, 6))
        grid = GridSpec(2, 5, figure=fig)

        for i, (correlation, position) in enumerate(most_correlated):
            section1 = most_corr_vectors[0][position:position + section_size]
            section2 = most_corr_vectors[1][position:position + section_size]

            plt.subplot(grid[0, i])
            plt.plot(section1, label='Hypermatrix Factor')
            plt.plot(section2, label='Insulation Score')
            plt.title(f'Most Correlated {i + 1}')
            plt.legend()

        for i, (correlation, position) in enumerate(least_correlated):
            section1 = most_corr_vectors[0][position:position + section_size]
            section2 = most_corr_vectors[1][position:position + section_size]

            plt.subplot(grid[1, i])
            plt.plot(section1, label='Hypermatrix Factor')
            plt.plot(section2, label='Insulation Score')
            plt.title(f'Least Correlated {i + 1}')
            plt.legend()

        # Adjust layout and save the figure to a file
        plt.tight_layout()
        plt.savefig(path + f'Images/res_{resolution}_window{insulation_window}_correlation_sections_ch{ch}_plot.png')
        # Close the figure 
        plt.close(fig)
correlatedSections
comment
