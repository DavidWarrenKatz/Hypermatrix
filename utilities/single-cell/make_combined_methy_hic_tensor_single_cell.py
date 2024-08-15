import h5py
import numpy as np
import os
from scipy.sparse import csr_matrix, identity
from sklearn.decomposition import NMF
from config_and_print import methy_directory, filtered_list, chrom_file, resolutions, output_directory
from scipy.ndimage import binary_dilation

def parse_resolution(resolution_str):
    if ':' in resolution_str:
        resolution_value, resolution_label = resolution_str.split(':')
        try:
            resolution = int(resolution_value)
            return resolution, resolution_label
        except ValueError:
            raise ValueError(f"Resolution value should be an integer: '{resolution_value}' in '{resolution_str}'")
    else:
        raise ValueError(f"Invalid resolution format: '{resolution_str}'. Expected format 'value:label', e.g., '1000000:1Mb'.")

def normalize_matrix_columns(A):
    column_norms = np.linalg.norm(A, axis=0)
    zero_norms = column_norms == 0
    column_norms[zero_norms] = 1
    normalized_A = A / column_norms
    normalized_A[:, zero_norms] = 0
    return normalized_A

def enforce_neighborhood_zeroing(methy_matrix, hic_matrix, neighborhood_size=4):
    # Create a boolean mask where hic_matrix is non-zero
    non_zero_mask = (hic_matrix != 0)
    
    # Dilate the non_zero_mask to include the neighborhood of the non-zero values
    dilated_non_zero_mask = binary_dilation(non_zero_mask, structure=np.ones((neighborhood_size*2 + 1, neighborhood_size*2 + 1)))
    
    # Set values in methy_matrix to zero wherever the dilated_non_zero_mask is False
    methy_matrix[~dilated_non_zero_mask] = 0
    
    return methy_matrix

def load_csr_matrix_from_hdf5(file_path):
    with h5py.File(file_path, 'r') as file:
        data = file['Matrix/data'][:]
        indices = file['Matrix/indices'][:]
        indptr = file['Matrix/indptr'][:]
        shape = file['Matrix'].attrs['shape']
    return csr_matrix((data, indices, indptr), shape=shape)

def normalize_matrix(A):
    norm = np.linalg.norm(A)
    if norm == 0:
        return A
    return A / norm

def csr_pearson_correlation(csr_mat):
    """Calculates Pearson correlation from a CSR matrix."""
    csc_mat = csr_mat.tocsc()
    mean = np.array(csc_mat.mean(axis=1)).flatten()
    std_dev = np.sqrt(csc_mat.power(2).mean(axis=1).A1 - mean**2)
    valid_std_dev = std_dev != 0
    rows, cols = csr_mat.nonzero()
    standardized_data = np.divide(csr_mat.data - mean[rows], std_dev[rows], where=valid_std_dev[rows])
    standardized_csr = csr_matrix((standardized_data, (rows, cols)), shape=csr_mat.shape)
    correlation_matrix = standardized_csr.dot(standardized_csr.T).toarray()
    diag = np.sqrt(np.diag(correlation_matrix))

    # Adjusting diagonal for valid standard deviations
    diag = np.where(valid_std_dev, diag, 1)  # replace zero with one to avoid division by zero
    correlation_matrix /= diag[:, None]
    correlation_matrix /= diag[None, :]
    return csr_matrix(np.nan_to_num(correlation_matrix))  # Replace NaNs with zero, caused by division by zero

def make_matrix_non_negative(matrix):
    # Convert the input to a NumPy array for easy manipulation
    matrix = np.array(matrix)
    
    # Find the minimum value in the matrix
    min_value = np.min(matrix)
    
    # Calculate the value needed to add to make all elements non-negative
    if min_value < 0:
        offset = -min_value
    else:
        offset = 0
    
    # Add the offset to every element in the matrix
    new_matrix = matrix + offset
    
    return new_matrix

def nmf(matrix, rank):
    model = NMF(n_components=rank, init='random', random_state=0)
    W = model.fit_transform(matrix)
    H = model.components_
    reconstructed_matrix = np.dot(W, H)
    return reconstructed_matrix

def ensure_same_size(matrix1, matrix2):
    """Ensure both matrices are the same size by padding with zeros as necessary."""
    max_rows = max(matrix1.shape[0], matrix2.shape[0])
    max_cols = max(matrix1.shape[1], matrix2.shape[1])

    def pad_matrix(matrix, new_shape):
        padded_matrix = np.zeros(new_shape)
        padded_matrix[:matrix.shape[0], :matrix.shape[1]] = matrix
        return padded_matrix

    matrix1 = pad_matrix(matrix1, (max_rows, max_cols))
    matrix2 = pad_matrix(matrix2, (max_rows, max_cols))

    return matrix1, matrix2

def process_tensor(methy_file_path, hic_path, output_tensor_path):
    try:
        with h5py.File(methy_file_path, 'r') as hf:
            dataset_name = list(hf.keys())[0]  
            methy_matrix = hf[dataset_name][:]
            print(f"Loaded methylation matrix with shape: {methy_matrix.shape}")
    
        # Check if Hi-C file exists and load the Hi-C matrix, otherwise create an identity matrix
        if os.path.exists(hic_path):
            hic_matrix = load_csr_matrix_from_hdf5(hic_path)
            hic_matrix_dense = hic_matrix.toarray()
            if hic_matrix_dense.shape[0] > 1 and hic_matrix_dense.shape[1] > 1:
                hic_matrix = hic_matrix_dense[:-1, :-1]
            else:
                hic_matrix = hic_matrix_dense
            print(f"Loaded Hi-C matrix with shape: {hic_matrix.shape}")
            
            # Ensure both matrices are the same size before setting elements to zero
            methy_matrix, hic_matrix = ensure_same_size(methy_matrix, hic_matrix)
            methy_matrix = enforce_neighborhood_zeroing(methy_matrix, hic_matrix, neighborhood_size=4)
            #methy_matrix[hic_matrix == 0] = 0
        else:
            hic_matrix = identity(methy_matrix.shape[0]).toarray()  # Convert to dense format
            print(f"Hi-C matrix not found, created identity matrix with shape: {hic_matrix.shape}")
            # Randomly zero out 90% of the methy_matrix entries
            random_mask = np.random.choice([0, 1], size=methy_matrix.shape, p=[0.9, 0.1])
            methy_matrix *= random_mask
            print(f"Applied random mask to methylation matrix")

        # Ensure both matrices are the same size
        methy_matrix, hic_matrix = ensure_same_size(methy_matrix, hic_matrix)
        print(f"Ensured matrices are the same size: {methy_matrix.shape} and {hic_matrix.shape}")
    
        methy_matrix = csr_matrix(methy_matrix)
        methy_matrix = make_matrix_non_negative(csr_pearson_correlation(methy_matrix).toarray())
        hic_matrix = csr_matrix(hic_matrix)
        hic_matrix = make_matrix_non_negative(csr_pearson_correlation(hic_matrix).toarray())
    
        # Normalize each matrix so that each column has a norm of 1
        # not sure if this is adding bias
        hic_matrix = normalize_matrix_columns(hic_matrix)
        methy_matrix = normalize_matrix_columns(methy_matrix)
    
        # Additional normalization step
        hic_matrix = normalize_matrix(hic_matrix)
        methy_matrix = normalize_matrix(methy_matrix)
    
        # Combine both normalized matrices into a tensor
        tensor = np.stack((hic_matrix, methy_matrix), axis=-1)
    
        # Create the output directory if it does not exist
        os.makedirs(os.path.dirname(output_tensor_path), exist_ok=True)
    
        # Save the tensor to an HDF5 file
        with h5py.File(output_tensor_path, 'w') as hf:
            hf.create_dataset('Tensor', data=tensor)
        print(f"Tensor saved to {output_tensor_path}")
    
    except Exception as e:
        print(f"Error processing tensor for {methy_file_path} and {hic_path}: {e}")

def generate_file_paths(chromosome, prefix, resolution_label, output_directory):
    methy_file_path = f'{output_directory}/methy_{resolution_label}_outerproduct_dir/{chromosome}/{prefix}_outer_product.h5'
    hic_path = f'{output_directory}/hic_{resolution_label}_emphasized_dir/{chromosome}/{prefix}_{chromosome}.h5'
    output_tensor_basepath = f'{output_directory}/hic_methy_{resolution_label}_tensor_singlecell/{chromosome}/'
    output_tensor_path = output_tensor_basepath + f'{prefix}_{chromosome}.h5'
    return methy_file_path, hic_path, output_tensor_path

# Ensure resolutions is treated as a tuple or list of strings
if isinstance(resolutions, str):
    resolutions = (resolutions,)

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]
resolution, resolution_label = parse_resolution(resolution_str)

print(f"Extracted resolution: {resolution}, label: {resolution_label}")

prefix_file_path = filtered_list
# Read prefixes from the file
with open(prefix_file_path, 'r') as f:
    prefixes = [line.strip() for line in f]    

for i in range(1, 23):  
    chromosome = f'chr{i}'
    for prefix in prefixes:
        methy_file_path, hic_path, output_tensor_path = generate_file_paths(chromosome, prefix, resolution_label, output_directory)
        
        if os.path.exists(output_tensor_path):
            print(f"File {output_tensor_path} already exists. Skipping.")
        else:
            # Ensure the output directory exists
            os.makedirs(os.path.dirname(output_tensor_path), exist_ok=True)
            process_tensor(methy_file_path, hic_path, output_tensor_path)
