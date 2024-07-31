import os
import numpy as np
import h5py
import pyBigWig
import scipy.stats
import matplotlib.pyplot as plt
import math
from config_and_print import data_path, resolutions, chromosomes, iterations

# Define paths for saved data
tensor_active_factor_path = data_path + 'Workspaces/genome_wide_tensor_active_factor.npy'
aiden_active_factor_path = data_path + 'Workspaces/genome_wide_aiden_active_factor.npy'
bin_averages_filtered_path = data_path + 'Workspaces/bin_averages_H3K4me3_filtered.npy'

def load_dark_bins(file_path):
    try:
        with h5py.File(file_path, 'r') as hf:
            dark_bins = hf['dark_bins_indices'][:]
        return dark_bins
    except Exception as e:
        print(f"[ERROR] Failed to load dark bins from {file_path}: {e}")
        return None

def remove_dark_bins(matrix, dark_bins):
    try:
        mask = np.ones(matrix.shape[0], dtype=bool)
        mask[dark_bins] = False
        filtered_matrix = matrix[mask][:, mask]
        return filtered_matrix
    except Exception as e:
        print(f"[ERROR] Failed to remove dark bins: {e}")
        return None
    
def normalize_vector(active_factor, threshold=1.5):
    mean = np.mean(active_factor)
    std = np.std(active_factor)
    z_scores = (active_factor - mean) / std
    active_factor = np.where(np.abs(z_scores) > threshold, mean + threshold * std * np.sign(z_scores), active_factor)
    norm = np.linalg.norm(active_factor)
    if norm == 0:
        raise ValueError("Cannot normalize a zero vector")
    normalized_vector = active_factor / norm
    return normalized_vector

def load_and_replace_nan(file_path):
    data = np.loadtxt(file_path)
    return np.nan_to_num(data)

def calculate_bin_averages(data, elements_per_bin):
    num_bins = math.ceil(len(data) / elements_per_bin)
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

def save_data_if_not_exists(file_path, data):
    if not os.path.exists(file_path):
        np.save(file_path, data)

if os.path.exists(tensor_active_factor_path) and os.path.exists(aiden_active_factor_path) and os.path.exists(bin_averages_filtered_path):
    genome_wide_tensor_active_factor = np.load(tensor_active_factor_path)
    genome_wide_aiden_active_factor = np.load(aiden_active_factor_path)
    genome_wide_H3K4me3 = np.load(bin_averages_filtered_path)
else:
    genome_wide_tensor_active_factor = []
    genome_wide_aiden_active_factor = []
    genome_wide_H3K4me3 = []

    for ch in chromosomes:
        for resolution in resolutions:
            tensor_factor_file = data_path +  f"Workspaces/individual/ch{ch}_res{resolution}_structedData_3rdCumulant_rank2_{iterations}iterations.h5"
            dataset_name = '/U'
            with h5py.File(tensor_factor_file, 'r') as f:
                U = f[dataset_name][:]
            tensor_factor1 = U[0, :]
            tensor_factor2 = U[1, :]
            tensor_factor3 = U[2, :]

            aiden_path = new_path + f'eigenvector/res{resolution}_ch{ch}_oe_KR_eigenvector.txt'
            aiden_eigenvector = load_and_replace_nan(aiden_path)

            dark_bins_file = new_path + f'Workspaces/individual/ch{ch}_res{resolution}_darkBins.h5'
            dark_bins = load_dark_bins(dark_bins_file)

            mask = np.ones(aiden_eigenvector.shape[0], dtype=bool)
            mask[dark_bins] = False
            filtered_aiden_eigenvector = aiden_eigenvector[mask]
            
            H3K4me3_file = "/home/dwk681/workspace/hypermatrix_test/hypermatrix/projects/softwarefiles/H3K4me3.bigwig"
            bw_H3K4me3 = pyBigWig.open(H3K4me3_file)
            start = 0
            end = bw_H3K4me3.chroms()[f'chr{ch}']
            data_H3K4me3 = bw_H3K4me3.values(f'chr{ch}', start, end)
            bin_averages = calculate_bin_averages(data_H3K4me3, resolution)
            mask = np.ones(bin_averages.shape[0], dtype=bool)
            mask[dark_bins] = False
            bin_averages_filtered = bin_averages[mask]

            correlation1 = np.corrcoef(tensor_factor1, bin_averages_filtered)[0, 1]
            correlation2 = np.corrcoef(tensor_factor2, bin_averages_filtered)[0, 1]
            correlation3 = np.corrcoef(tensor_factor3, bin_averages_filtered)[0, 1]
            if correlation1 > correlation2 and correlation1 > correlation3:
                active_factor = tensor_factor1
            elif correlation2 > correlation1 and correlation2 > correlation3:
                active_factor = tensor_factor2
            else:
                active_factor = tensor_factor3

            correlation1 = np.corrcoef(filtered_aiden_eigenvector, bin_averages_filtered)[0, 1]
            correlation2 = np.corrcoef(-filtered_aiden_eigenvector, bin_averages_filtered)[0, 1]
            if correlation1 > correlation2:
                filtered_aiden_eigenvector = filtered_aiden_eigenvector
            else:
                filtered_aiden_eigenvector = -filtered_aiden_eigenvector

            active_factor = normalize_vector(active_factor)

            genome_wide_tensor_active_factor.extend(active_factor)
            genome_wide_aiden_active_factor.extend(filtered_aiden_eigenvector)
    
    genome_wide_tensor_active_factor = np.array(genome_wide_tensor_active_factor)
    genome_wide_aiden_active_factor = np.array(genome_wide_aiden_active_factor)
    save_data_if_not_exists(tensor_active_factor_path, genome_wide_tensor_active_factor)
    save_data_if_not_exists(aiden_active_factor_path, genome_wide_aiden_active_factor)
    save_data_if_not_exists(bin_averages_filtered_path, bin_averages_filtered)
    
# Create a density scatter plot
x = genome_wide_tensor_active_factor
y = genome_wide_aiden_active_factor
xy = np.vstack([x, y])
z = scipy.stats.gaussian_kde(xy)(xy)

fig, ax = plt.subplots(figsize=(15, 7))
sc = ax.scatter(x, y, c=z, s=10, cmap='viridis')
fig.colorbar(sc, label='Density')
ax.set_title(f'3rd Cumulant Hi-C Tensor A Compartment Calls\nRes {resolution}', fontsize=20)
ax.set_xlabel('Tensor Active Factor')
ax.set_ylabel('Aiden Eigenvector')

correlation_coefficient_genome_wide, p_value_genome_wide = scipy.stats.pearsonr(x, y)
caption = f'The first two factors of a rank two decomposition are always negatively correlated.\nCh1-22 Average Correlation: {correlation_coefficient_genome_wide:.2f}, average p-value: {p_value_genome_wide:.2f}.'
fig.text(0.5, 0.02, caption, ha='center', fontsize=12)

plt.tight_layout(rect=[0, 0.1, 1, 0.95])
plt.savefig('density_scatter_plot.png')
plt.show()

