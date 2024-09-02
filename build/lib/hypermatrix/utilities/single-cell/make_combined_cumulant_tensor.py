import h5py
import numpy as np
import os
from config_and_print import filtered_list, chrom_file, resolutions, output_directory

# Ensure resolutions is treated as a tuple or list of strings
if isinstance(resolutions, str):
    resolutions = (resolutions,)

# Extract resolution value and label from the resolutions string
resolution_str = resolutions[0]

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

resolution, resolution_label = parse_resolution(resolution_str)
prefix_file_path = filtered_list

def load_cumulant_tensor(file_path):
    """Load the cumulant tensor from an HDF5 file."""
    with h5py.File(file_path, 'r') as h5file:
        cumulant_tensor = h5file['cumulant_tensor'][:]
    return cumulant_tensor

def make_tensor_non_negative(tensor):
    """Ensure that the tensor is non-negative by adding the smallest value if it's negative."""
    min_value = tensor.min()
    if min_value < 0:
        tensor += np.abs(min_value)
    return tensor

def normalize_cumulant_tensor(tensor, positive_threshold, negative_threshold):
    """Normalize the cumulant tensor based on the specified thresholds."""
    normalized_tensor = np.zeros_like(tensor)

    # Binarize positive values
    normalized_tensor[tensor > positive_threshold] = positive_threshold
    
    # Binarize negative values
    normalized_tensor[(tensor < 0) & (np.abs(tensor) > negative_threshold)] = -1

    return normalized_tensor

def crop_tensor_to_match(reference_tensor, tensor_to_crop):
    """Crop the tensor_to_crop to match the dimensions of the reference_tensor."""
    ref_shape = reference_tensor.shape
    crop_shape = tensor_to_crop.shape
    
    # Determine the minimum shape for each dimension
    min_shape = [min(ref_shape[i], crop_shape[i]) for i in range(3)]
    
    # Crop the tensors to the minimum shape
    cropped_tensor = tensor_to_crop[:min_shape[0], :min_shape[1], :min_shape[2]]
    
    return cropped_tensor

def adjust_methy_tensor_to_hic(hic_tensor, methy_tensor, neighborhood_size=1):
    """Adjust the methylation cumulant tensor to match the non-zero pattern of the Hi-C cumulant."""
    adjusted_methy_tensor = np.zeros_like(methy_tensor)

    nonzero_indices = np.argwhere(hic_tensor != 0)
    
    for idx in nonzero_indices:
        x, y, z = idx
        # Apply a neighborhood around the non-zero entries
        x_min, x_max = max(0, x-neighborhood_size), min(methy_tensor.shape[0], x+neighborhood_size+1)
        y_min, y_max = max(0, y-neighborhood_size), min(methy_tensor.shape[1], y+neighborhood_size+1)
        z_min, z_max = max(0, z-neighborhood_size), min(methy_tensor.shape[2], z+neighborhood_size+1)
        adjusted_methy_tensor[x_min:x_max, y_min:y_max, z_min:z_max] = methy_tensor[x_min:x_max, y_min:y_max, z_min:z_max]

    return adjusted_methy_tensor

def normalize_tensor(tensor):
    """Normalize a tensor to have a norm of one."""
    norm = np.linalg.norm(tensor)
    if norm > 0:
        return tensor / norm
    return tensor

def process_and_combine_tensors(chromosome, prefix_hic, prefix_methy, positive_threshold, negative_threshold, neighborhood_size, output_directory):
    hic_file_path = os.path.join(output_directory, f'hic_{resolution_label}_cumulant_dir', chromosome, f'{prefix_hic}_{chromosome}_cumulant.h5')
    methy_file_path = os.path.join(output_directory, f'methy_{resolution_label}_cumulant_dir', chromosome, f'{prefix_methy}_cumulant.h5')

    # Check if Hi-C cumulant file exists
    if not os.path.exists(hic_file_path):
        print(f"Hi-C cumulant file does not exist for {prefix_hic} on {chromosome}. Skipping.")
        return

    # Check if methylation cumulant file exists
    if not os.path.exists(methy_file_path):
        print(f"Methylation cumulant file does not exist for {prefix_methy} on {chromosome}. Skipping.")
        return

    # Load tensors
    hic_tensor = load_cumulant_tensor(hic_file_path)
    methy_tensor = load_cumulant_tensor(methy_file_path)

    # Ensure both tensors have the same dimensions by cropping the larger tensor
    if hic_tensor.shape != methy_tensor.shape:
        if hic_tensor.shape > methy_tensor.shape:
            hic_tensor = crop_tensor_to_match(methy_tensor, hic_tensor)
        else:
            methy_tensor = crop_tensor_to_match(hic_tensor, methy_tensor)

    # Normalize the Hi-C cumulant tensor
    hic_tensor = normalize_cumulant_tensor(hic_tensor, positive_threshold, negative_threshold)

    # Adjust the methylation tensor to match the Hi-C tensor's non-zero pattern
    methy_tensor = adjust_methy_tensor_to_hic(hic_tensor, methy_tensor, neighborhood_size)

    # Normalize both tensors to have a norm of one
    hic_tensor = normalize_tensor(hic_tensor)
    methy_tensor = normalize_tensor(methy_tensor)

    # Ensure both tensors are non-negative
    hic_tensor = make_tensor_non_negative(hic_tensor)
    methy_tensor = make_tensor_non_negative(methy_tensor)

    # Combine the tensors into a single tensor
    combined_tensor = np.stack((hic_tensor, methy_tensor), axis=-1)

    # Prepare the output directory
    combined_output_dir = os.path.join(output_directory, f'{resolution_label}_combined_cumulant', chromosome)
    os.makedirs(combined_output_dir, exist_ok=True)

    # Save the combined tensor
    combined_file_path = os.path.join(combined_output_dir, f'{prefix_hic}_{chromosome}_combined_cumulant.h5')
    with h5py.File(combined_file_path, 'w') as h5file:
        h5file.create_dataset('combined_cumulant_tensor', data=combined_tensor)
    print(f"Combined cumulant tensor saved to {combined_file_path}")

chromosomes = [f'chr{i}' for i in range(1, 23)]  # Chromosomes 1 to 22

# Read prefixes from the file
with open(prefix_file_path, 'r') as f:
    prefixes = [line.strip() for line in f]

positive_threshold = 1  # Adjust as needed
negative_threshold = 0.5  # Adjust as needed
neighborhood_size = 1  # Size of neighborhood around non-zero Hi-C entries

# Process and combine tensors for each chromosome and prefix
for chromosome in chromosomes:
    for prefix in prefixes:
        process_and_combine_tensors(chromosome, prefix, prefix, positive_threshold, negative_threshold, neighborhood_size, output_directory)

