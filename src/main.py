import subprocess
import sys
import os

# resolving filepaths locally or in project folder
# upgrading parsing arguments
# create submodules



def check_and_install_requirements():
    try:
        import tensorly
    except ImportError:
        print("Required packages not found. Installing...")
        # Get the directory of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Construct the path to installation.py
        installation_script = os.path.join(current_dir, 'installation.py')
        # Check if the installation script exists
        if not os.path.isfile(installation_script):
            sys.exit("Error: installation.py script not found.")
        subprocess.check_call([sys.executable, installation_script])

def normalize_vector(active_factor):
    import numpy as np
    norm = np.linalg.norm(active_factor)
    if norm == 0:
        raise ValueError("Cannot normalize a zero vector")
    return active_factor / norm

def load_and_replace_nan(file_path):
    import numpy as np
    data = np.loadtxt(file_path)
    return np.nan_to_num(data, nan=0.0)

def calculate_bin_averages(data, elements_per_bin):
    import numpy as np
    import math
    num_bins = math.ceil(len(data) / elements_per_bin)
    bin_averages = np.zeros(num_bins)
    for i in range(num_bins):
        start_index = i * elements_per_bin
        end_index = min((i + 1) * elements_per_bin, len(data))
        bin_data = data[start_index:end_index]
        bin_averages[i] = np.mean(bin_data) if len(bin_data) > 0 else 0
    return bin_averages

def process_data(path, resolution):
    import scipy.io
    import os
    import numpy as np
    import pyBigWig

    sol_factors_downsampled_by_ch = {}
    components_by_channel = {ch: [] for ch in range(1, 23)}
    removed_rowIndices_by_ch = {ch: [] for ch in range(1, 23)}
    filtered_factor1 = []
    filtered_factor2 = []
    filtered_factor3 = []

    for ch in range(1, 23):
        mat_file = os.path.join(path, f'Workspaces/individual/ch{ch}_res{resolution}_structedData_10Downsampled_pearsons_first3_500iterations.mat')
        data = scipy.io.loadmat(mat_file)
        sol_factors_downsampled_by_ch[ch] = data['sol_factors_downsampled']

        removed_rows_path = os.path.join(path, f'Workspaces/individual/ch{ch}_res{resolution}_observed_NONE_removedRows.mat')
        data2 = scipy.io.loadmat(removed_rows_path)
        removed_rowIndices_by_ch[ch] = data2['removed_rows_indices']

        N_ch = sol_factors_downsampled_by_ch[ch][0, 1][0, 0][:, 0].size
        filtered_indices_ch = [i for i in range(N_ch-2) if i not in removed_rowIndices_by_ch[ch]]

        for component in range(1, 3):
            component_data = data['sol_factors_downsampled'][0, 1][0, 0][:, component - 1]
            components_by_channel[ch].append(component_data)

        vector1 = components_by_channel[ch][0][filtered_indices_ch]
        vector2 = components_by_channel[ch][1][filtered_indices_ch]

        bigwigfile1 = os.path.join(path, 'ENCFF154XCY_hg19_H3K4me3_GM12878.bigWig')
        bw_file1 = pyBigWig.open(bigwigfile1)
        start = 0
        end = bw_file1.chroms()[f'chr{ch}']
        data1 = bw_file1.values(f'chr{ch}', start, end)
        bin_averages1 = calculate_bin_averages(data1, resolution)
        bin_averages1_filtered = bin_averages1[filtered_indices_ch]

        correlation1 = np.corrcoef(vector1, bin_averages1_filtered)[0, 1]
        correlation2 = np.corrcoef(vector2, bin_averages1_filtered)[0, 1]
        active_factor = vector1 if correlation1 > correlation2 else vector2
        passive_factor = vector2 if correlation1 > correlation2 else vector1

        file_path = os.path.join(path, f'eigenvector/res{resolution}_ch{ch}_observed_NONE_NONE_eigenvector.txt')
        vector3 = load_and_replace_nan(file_path)[filtered_indices_ch]
        vector4 = -vector3
        correlation1 = np.corrcoef(vector3, bin_averages1_filtered)[0, 1]
        correlation2 = np.corrcoef(vector4, bin_averages1_filtered)[0, 1]
        filtered_aiden_eigenvectors = vector3 if correlation1 > correlation2 else vector4

        active_factor = normalize_vector(active_factor)
        passive_factor = normalize_vector(passive_factor)

        filtered_factor1.extend(active_factor)
        filtered_factor2.extend(passive_factor)
        filtered_factor3.extend(filtered_aiden_eigenvectors)

    return filtered_factor1, filtered_factor2, filtered_factor3, components_by_channel

def plot_results(filtered_factor1, filtered_factor2, filtered_factor3, components_by_channel, chromosome, resolution, path):
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    fig.suptitle(f'Factors from Downsampled Hi-C Data\nres {resolution}', fontsize=20)

    filtered_indices_ch = [i for i in range(len(components_by_channel[chromosome][0]))]

    correlation_coefficient, p_value = pearsonr(components_by_channel[chromosome][0][filtered_indices_ch], components_by_channel[chromosome][1][filtered_indices_ch])
    file_path = os.path.join(path, f'eigenvector/res{resolution}_ch{chromosome}_observed_NONE_NONE_eigenvector.txt')
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

    vmin, vmax = np.min(filtered_factor3), np.max(filtered_factor3)
    correlation_coefficient_genome_wide, p_value_genome_wide = pearsonr(filtered_factor1, filtered_factor2)
    sc = ax2.scatter(filtered_factor1, filtered_factor2, c=filtered_factor3, cmap='viridis', vmin=vmin, vmax=vmax)
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
    correlation_coefficient, p_value = pearsonr(components_by_channel[chromosome][0][filtered_indices_ch], components_by_channel[chromosome][1][filtered_indices_ch])
    ax.plot(components_by_channel[chromosome][0], 'b', linewidth=2, label='Factor One')
    ax.plot(components_by_channel[chromosome][1], 'r', linewidth=2, label='Factor Two')
    ax.set_xlabel('Genomic Bin')
    ax.set_ylabel('Value')
    ax.set_title(f'Plot of chr {chromosome}\nCorrelation:{correlation_coefficient:.2f}, p-value: {p_value:.2f}')
    ax.legend()
    ax.grid(True)

    plt.show()

def process_and_plot(path, resolution, chromosome):
    filtered_factor1, filtered_factor2, filtered_factor3, components_by_channel = process_data(path, resolution)
    plot_results(filtered_factor1, filtered_factor2, filtered_factor3, components_by_channel, chromosome, resolution, path)

if __name__ == '__main__':
    check_and_install_requirements()

    import argparse

    parser = argparse.ArgumentParser(description='Process Hi-C data and visualize the results.')
    parser.add_argument('--path', type=str, required=True, help='Path to the data files')
    parser.add_argument('--resolution', type=int, default=1_000_000, help='Resolution of the data')
    parser.add_argument('--chromosome', type=int, default=11, help='Chromosome to visualize')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    process_and_plot(args.path, args.resolution, args.chromosome)
