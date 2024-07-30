import numpy as np
import h5py
import os
import sys
from sklearn.decomposition import NMF


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

def compute_correlation_matrix(matrix):
    try:
        correlation_matrix = np.corrcoef(matrix)
        min_value = np.min(np.nan_to_num(correlation_matrix))
        if min_value < 0:
            correlation_matrix += np.abs(min_value)

        return correlation_matrix
    except Exception as e:
        print(f"[ERROR] Failed to compute correlation matrix: {e}")
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
                    normalized_kappa_ijk[i, j, k] = kappa_ijk[i, j, k] / denominator if denominator != 0 else 0

        return normalized_kappa_ijk
    except Exception as e:
        print(f"[ERROR] Failed to compute third-order cumulant: {e}")
        return None


def compute_degree_2_cumulant(matrix):
    try:
        n, m = matrix.shape
        mean_vec = np.mean(matrix, axis=1)
        std_vec = np.std(matrix, axis=1)

        degree_2_cumulant = np.zeros((n, n, n))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    # Calculate numerator
                    numerator = np.sum((matrix[i] - mean_vec[i]) * (matrix[j] - mean_vec[j]) * (matrix[k] - mean_vec[k]))
                    
                    # Calculate denominator
                    denom_i = np.sum((matrix[i] - mean_vec[i])**2)
                    denom_j = np.sum((matrix[j] - mean_vec[j])**2)
                    denom_k = np.sum((matrix[k] - mean_vec[k])**2)
                    denominator = np.sqrt(denom_i * denom_j * denom_k)
                    
                    # Compute cumulant with NaN handling
                    if denominator != 0:
                        degree_2_cumulant[i, j, k] = numerator / denominator
                    else:
                        degree_2_cumulant[i, j, k] = np.nan

        min_value = np.min(np.nan_to_num(degree_2_cumulant))
        if min_value < 0:
            degree_2_cumulant -= min_value


        return degree_2_cumulant
    except Exception as e:
        print(f"[ERROR] Failed to compute degree-2 cumulant: {e}")
        return None

def compute_non_negative_rank_2_decomposition(matrix):
    try:
        model = NMF(n_components=2, init='random', random_state=0)
        W = model.fit_transform(matrix)
        H = model.components_
        return W, H
    except Exception as e:
        print(f"[ERROR] Failed to compute non-negative rank-2 decomposition: {e}")
        return None, None

def process_hic_files(path, resolutions, chromosomes, data_types):
    for resolution in resolutions:
        for chromosome in chromosomes:
            for data_type in data_types:
                hic_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR.h5'
                dark_bins_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_darkBins.h5'

                hic_matrix = load_hic_matrix(hic_file)
                dark_bins = load_dark_bins(dark_bins_file)

                if hic_matrix is not None and dark_bins is not None:
                    filtered_matrix = remove_dark_bins(hic_matrix, dark_bins)

                    if filtered_matrix is not None:
                        corr_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_corr_new.h5'
                        cumulant_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_cumulant_new.h5'
                        nn_decomp_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_nn_decomp.h5'

                        # Check if correlation matrix file already exists
                        if not os.path.exists(corr_output_file):
                            correlation_matrix = compute_correlation_matrix(filtered_matrix)
                            if correlation_matrix is not None:
                                try:
                                    with h5py.File(corr_output_file, 'w') as hf:
                                        hf.create_dataset('correlation_matrix', data=correlation_matrix)
                                    print(f'[INFO] Saved correlation matrix to {corr_output_file}')
                                except Exception as e:
                                    print(f"[ERROR] Failed to save correlation matrix: {e}")
                        else:
                            print(f'[INFO] Correlation matrix already exists: {corr_output_file}')

                        # Check if degree-2 cumulant matrix file already exists
                        if not os.path.exists(cumulant_output_file):
                            degree_2_cumulant = compute_third_order_cumulant(filtered_matrix)
                            if degree_2_cumulant is not None:
                                try:
                                    with h5py.File(cumulant_output_file, 'w') as hf:
                                        hf.create_dataset('degree_2_cumulant', data=degree_2_cumulant)
                                    print(f'[INFO] Saved degree-2 cumulant matrix to {cumulant_output_file}')
                                except Exception as e:
                                    print(f"[ERROR] Failed to save degree-2 cumulant matrix: {e}")
                        else:
                            print(f'[INFO] Degree-2 cumulant matrix already exists: {cumulant_output_file}')

                        if not os.path.exists(nn_decomp_output_file):
                            W, H = compute_non_negative_rank_2_decomposition(correlation_matrix)
                            if W is not None and H is not None:
                                save_matrix_to_file(W, nn_decomp_output_file, 'W')
                                save_matrix_to_file(H, nn_decomp_output_file, 'H')



if __name__ == "__main__":
    data_path = sys.argv[1]    
    resolutions_list = [int(r.strip().replace("'", "")) for r in sys.argv[2].split(",")]
    chromosomes_list = [c.strip().replace("'", "") for c in sys.argv[3].split(",")]
    data_types_list = [d.strip().replace("'", "") for d in sys.argv[4].split(",")]
    process_hic_files(data_path, resolutions_list, chromosomes_list, data_types_list)
