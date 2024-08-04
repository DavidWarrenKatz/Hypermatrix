import numpy as np
import h5py
import os
import sys
from config_and_print import mappability_threshold

def load_hic_matrix(file_path):
    try:
        with h5py.File(file_path, 'r') as hf:
            matrix = hf['matrix'][:]
        return matrix
    except Exception as e:
        print(f"[ERROR] Failed to load Hi-C matrix from {file_path}: {e}")
        return None

def load_dark_bins(file_path):
    try:
        with h5py.File(file_path, 'r') as hf:
            dark_bins = hf['dark_bins_indices'][:]
        return dark_bins
    except Exception as e:
        print(f"[ERROR] Failed to load dark bins from {file_path}: {e}")
        return None

def save_matrix_to_file(matrix, file_path, dataset_name):
    try:
        with h5py.File(file_path, 'w') as hf:
            hf.create_dataset(dataset_name, data=matrix)
        print(f'[INFO] Saved {dataset_name} to {file_path}')
    except Exception as e:
        print(f"[ERROR] Failed to save {dataset_name}: {e}")

def remove_dark_bins(matrix, dark_bins):
    try:
        mask = np.ones(matrix.shape[0], dtype=bool)
        mask[dark_bins] = False
        filtered_matrix = matrix[mask][:, mask]
        return filtered_matrix
    except Exception as e:
        print(f"[ERROR] Failed to remove dark bins: {e}")
        return None

# Define the function to compute the third-order cumulant
def compute_third_order_cumulant(matrix):
    try:
        n, m = matrix.shape
        mean_vec = np.mean(matrix, axis=1, keepdims=True)
        centered_matrix = matrix - mean_vec
        m_ijk = np.zeros((n, n, n))

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    m_ijk[i, j, k] = np.mean(centered_matrix[i] * centered_matrix[j] * centered_matrix[k])
        
        m_ij = np.dot(centered_matrix, centered_matrix.T) / m
        m_i = np.mean(centered_matrix, axis=1)

        kappa_ijk = np.zeros((n, n, n))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    kappa_ijk[i, j, k] = (
                        m_ijk[i, j, k]
                        - (m_i[i] * m_ij[j, k] + m_i[j] * m_ij[i, k] + m_i[k] * m_ij[i, j])
                        + 2 * m_i[i] * m_i[j] * m_i[k]
                    )

        # Normalize by dividing by the standard deviations
        std_vec = np.std(matrix, axis=1, keepdims=True)
        normalized_kappa_ijk = np.zeros((n, n, n))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    denominator = std_vec[i] * std_vec[j] * std_vec[k]
                    normalized_kappa_ijk[i, j, k] = kappa_ijk[i, j, k] / denominator if denominator != 0 else np.nan

        min_value = np.min(np.nan_to_num(normalized_kappa_ijk))
        if min_value < 0:
            normalized_kappa_ijk -= min_value

        return normalized_kappa_ijk
    except Exception as e:
        print(f"[ERROR] Failed to compute third-order cumulant: {e}")
        return None

def process_hic_files(path, resolutions, chromosomes, data_types):
    for resolution in resolutions:
        for chromosome in chromosomes:
            for data_type in data_types:
                hic_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR.h5'
                dark_bins_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_darkBins_mappability{mappability_threshold}.h5'

                hic_matrix = load_hic_matrix(hic_file)
                dark_bins = load_dark_bins(dark_bins_file)

                if hic_matrix is not None and dark_bins is not None:
                    filtered_matrix = remove_dark_bins(hic_matrix, dark_bins)

                    if filtered_matrix is not None:
                        cumulant_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_cumulant.h5'
                        # Check if degree-2 cumulant matrix file already exists
                        if not os.path.exists(cumulant_output_file):
                            degree_3_cumulant = compute_third_order_cumulant(filtered_matrix)
                            if degree_3_cumulant is not None:
                                try:
                                    with h5py.File(cumulant_output_file, 'w') as hf:
                                        hf.create_dataset('degree_3_cumulant', data=degree_3_cumulant)
                                    print(f'[INFO] Saved degree-3 cumulant matrix to {cumulant_output_file}')
                                except Exception as e:
                                    print(f"[ERROR] Failed to save degree-3 cumulant matrix: {e}")
                        else:
                            print(f'[INFO] Degree-3 cumulant matrix already exists: {cumulant_output_file}')

if __name__ == "__main__":
    data_path = sys.argv[1]    
    resolutions_list = [int(r.strip().replace("'", "")) for r in sys.argv[2].split(",")]
    chromosomes_list = [c.strip().replace("'", "") for c in sys.argv[3].split(",")]
    data_types_list = [d.strip().replace("'", "") for d in sys.argv[4].split(",")]
    process_hic_files(data_path, resolutions_list, chromosomes_list, data_types_list)
