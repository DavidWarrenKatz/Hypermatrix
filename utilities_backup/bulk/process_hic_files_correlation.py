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

def compute_non_negative_rank_2_decomposition(matrix):
    try:
        # Replace NaNs with zeros
        matrix = np.nan_to_num(matrix, nan=0.0)
        
        # Apply non-negative matrix factorization
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
                        corr_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_corr.h5'
                        cumulant_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_cumulant.h5'
                        #nn_decomp_output_file = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data_type}_KR_nn_decomp.h5'

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

                        if not os.path.exists(nn_decomp_output_file):
                            correlation_matrix = compute_correlation_matrix(filtered_matrix)
                            W, H = compute_non_negative_rank_2_decomposition(correlation_matrix)
                            if W is not None:
                                save_matrix_to_file(W, nn_decomp_output_file, 'W')
                                

if __name__ == "__main__":
    data_path = sys.argv[1]    
    resolutions_list = [int(r.strip().replace("'", "")) for r in sys.argv[2].split(",")]
    chromosomes_list = [c.strip().replace("'", "") for c in sys.argv[3].split(",")]
    data_types_list = [d.strip().replace("'", "") for d in sys.argv[4].split(",")]
    process_hic_files(data_path, resolutions_list, chromosomes_list, data_types_list)
