import os
import numpy as np
import h5py
import pyBigWig
from config_and_print import data_path, resolutions, chromosomes

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

def normalize_vector(active_factor, threshold=3):
    # Calculate the z-score for each element in the vector
    mean = np.mean(active_factor)
    std = np.std(active_factor)
    z_scores = (active_factor - mean) / std

    # Threshold the outliers
    active_factor = np.where(np.abs(z_scores) > threshold, mean + threshold * std * np.sign(z_scores), active_factor)

    # Normalize the vector
    norm = np.linalg.norm(active_factor)
    if norm == 0:
        raise ValueError("Cannot normalize a zero vector")
    normalized_vector = active_factor / norm
    return normalized_vector

def load_and_replace_nan(file_path):
    data = np.loadtxt(file_path)
    return np.nan_to_num(data)

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

sol_factors_downsampled_by_ch = {}
components_by_channel = {ch: [] for ch in range(1, 23)}
removed_rowIndices_by_ch = {ch: [] for ch in range(1, 23)}
filtered_factor1 = []
filtered_factor2 = []
filtered_factor3 = []

for ch in chromosomes:
    for resolution in resolutions:
        #First, load the three tensor factors for this chromosome
        tensor_factor_file = data_path + f'Workspaces/individual/ch{ch}_res{resolution}_structedData_3rdCumulant_rank3_400iterations.h5'
        dataset_name = '/U'
        # Open the HDF5 file and read the data
        with h5py.File(file_path, 'r') as f:
            U = f[dataset_name][:]
        tensor_factor1 = U[0,:]
        tensor_factor2 = U[1,:]
        tensor_factor3 = U[2,:]    

        #Second, load the Aiden Eigenvector for this chromosome for comparison
        aiden_path = data_path + f'eigenvector/res{resolution}_ch{chromosome}_observed_NONE_NONE_eigenvector.txt'
        aiden_eigenvector = load_and_replace_nan(aiden_path)
            
        #Third, load the dark bins, which were already removed for the tensor factors
        dark_bins_file = data_path + f'Workspaces/individual/ch{chromosome}_res{resolution}_darkBins.h5'    
        dark_bins = load_dark_bins(dark_bins_file)

        #Fourth, remove dark bins from Aiden eigenvector
        mask = np.ones(aiden_eigenvector.shape[0], dtype=bool)
        mask[dark_bins] = False
        filtered_aiden_eigenvector = aiden_eigenvector[mask]

        #Fifth, laod H3K4 file and remove dark bins
        H3K4me3_file="../projects/softwarefiles/H3K4me3.bigwig"
        bw_H3K4me3 = pyBigWig.open(H3K4me3_file)
        start = 0
        end = bw_H3K4me3.chroms()[chromosome_string]
        data_H3K4me3 = bw_H3K4me3.values(ch, start, end)
        bin_averages = calculate_bin_averages(data_H3K4me3, resolution)
        mask = np.ones(bin_averages.shape[0], dtype=bool)
        mask[dark_bins] = False
        bin_averages_filtered = bin_averages[mask]

        #Sixth, determine which of the three factors most correlates with H3K4me3
        #This will be the A compartment factor
        correlation1 = np.corrcoef(tensor_factor1, bin_averages_filtered)[0, 1]
        correlation2 = np.corrcoef(tensor_factor2, bin_averages_filtered)[0, 1]
        correlation3 = np.corrcoef(tensor_factor3, bin_averages_filtered)[0, 1]
        # Determine which factor has the highest correlation
        if correlation1 > correlation2 and correlation1 > correlation3:
            active_factor = tensor_factor1
        elif correlation2 > correlation1 and correlation2 > correlation3:
            active_factor = tensor_factor2
        else:
            active_factor = tensor_factor3

        #Seventh, determine if filtered_aiden_eigenvector or its negative is the active factor
        correlation1 = np.corrcoef(filtered_aiden_eigenvector, bin_averages_filtered)[0, 1]
        correlation2 = np.corrcoef(-filtered_aiden_eigenvector, bin_averages_filtered)[0, 1]
        if correlation1 > correlation2:
            filtered_aiden_eigenvector = filtered_aiden_eigenvector
        else:
            filtered_aiden_eigenvector = -filtered_aiden_eigenvector

        #normalize so that the norm of the tensor active factor is 1
        #this accounts for the arbitryness of only being defines up to a multiplicative constant
        #also remove outliers
        active_factor = normalize_vector(active_factor)
        passive_factor = normalize_vector(passive_factor)

        filtered_factor1.extend(active_factor)
        filtered_factor2.extend(passive_factor)
        filtered_factor3.extend(filtered_aiden_eigenvectors)

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
fig.suptitle(f'Factors from Downsampled Hi-C Data\nres {resolution}', fontsize=20)

#axis1
#A and B factors colored by AidenLab's eigenvector
correlation_coefficient, p_value = scipy.stats.pearsonr(components_by_channel[chromosome][0][filtered_indices_ch], components_by_channel[chromosome][1][filtered_indices_ch])
file_path = path + f'eigenvector/res{resolution}_ch{chromosome}_observed_NONE_NONE_eigenvector.txt'
vector3 = load_and_replace_nan(file_path)[filtered_indices_ch]
vmin, vmax = np.min(vector3), np.max(vector3)
sc = ax1.scatter(components_by_channel[chromosome][0][filtered_indices_ch], components_by_channel[chromosome][1][filtered_indices_ch], c=vector3, cmap='viridis', vmin=vmin, vmax=vmax)
ax1.set_title(f'Scatter Plot of A and B compartment Factors\nchr {chromosome}\nCorrelation:{correlation_coefficient:.2f}, p-value: {p_value:.2f}\nColored By Aiden Eigenvector Value')
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax1)
cbar.set_label('Aiden Eigenvector')
ax1.set_xlabel('Factor 1')
ax1.set_ylabel('Factor 2')

#axis2
#genome_wide figure of A and B Factors colored by Aiden Lab's Eigenvector
vmin, vmax = np.min(filtered_factor3), np.max(filtered_factor3)
correlation_coefficient_genome_wide, p_value_genome_wide = scipy.stats.pearsonr(filtered_factor1, filtered_factor2)
sc = ax2.scatter(filtered_factor1, filtered_factor2, c=filtered_factor3, cmap='viridis')
ax2.set_title(f'Autosomal Chromosome Scatter Plot\nColored By Aiden Eigenvector Value\nCorrelation:{correlation_coefficient_genome_wide:.2f}, p-value: {p_value_genome_wide:.2f}')
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax2)
cbar.set_label('Aiden Eigenvector')
ax2.set_xlabel('Factor 1')
ax2.set_ylabel('Factor 2')

caption = f'The first two factors of a rank two decomposition are always negatively correlated.\nCh1-22 Average Correlation: {correlation_coefficient_genome_wide:.2f}, average p-value: {p_value_genome_wide:.2f}.'
fig.text(0.5, 0.02, caption, ha='center', fontsize=12)
plt.tight_layout(rect=[0, 0.1, 1, 0.95])

fig2, ax = plt.subplots()
correlation_coefficient, p_value = scipy.stats.pearsonr(components_by_channel[chromosome][0][filtered_indices_ch], components_by_channel[chromosome][1][filtered_indices_ch])
ax.plot(components_by_channel[chromosome][0], 'b', linewidth=2, label='Factor One')
ax.plot(components_by_channel[chromosome][1], 'r', linewidth=2,  label='Factor Two')
ax.set_xlabel('Genomic Bin')
ax.set_ylabel('Value')
ax.set_title(f'Plot of chr {chromosome}\nCorrelation:{correlation_coefficient:.2f}, p-value: {p_value:.2f}')
ax.legend()
ax.grid(True)

plt.show()

#rename the factors for later usage
bulk_hic_a_compartment_filtered = filtered_factor1 
bulk_hic_b_compartment_filtered = filtered_factor2
bulk_hic_aiden_compartment_filtered = filtered_factor3
