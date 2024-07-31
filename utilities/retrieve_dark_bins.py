import os
import numpy as np
import h5py
import pyBigWig
from config_and_print import data_path, resolutions, chromosomes, mappability_threshold

def get_dark_bins(chromosome, resolution, bigwig):
    chrom_name = f'chr{chromosome}'
    if chrom_name not in bigwig.chroms():
        raise ValueError(f"Chromosome {chrom_name} not found in BigWig file.")
    
    chrom_size = bigwig.chroms(chrom_name)
    if chrom_size is None:
        raise ValueError(f"Size for chromosome {chrom_name} is None.")
    
    num_bins = chrom_size // resolution + (1 if chrom_size % resolution else 0)
    dark_bins = []

    for i in range(num_bins):
        start = i * resolution
        end = min(start + resolution, chrom_size)  # Ensure end does not exceed chrom_size
        values = bigwig.values(chrom_name, start, end, numpy=True)
        if len(values) == 0:
            continue
        mean_value = np.nanmean(values)
        if mean_value < mappability_threshold:
            dark_bins.append(i)
    
    return dark_bins


# Open the BigWig file for dark regions
dark_regions_file = "../projects/softwarefiles/dark_regions_hg19.bigWig"
bw = pyBigWig.open(dark_regions_file)

for resolution in resolutions:
    for chromosome in chromosomes:
        # Construct the HDF5 filename
        h5_filename = data_path + f'Workspaces/individual/ch{chromosome}_res{resolution}_darkBins_mappability{mappability_threshold}.h5'
        
        # Check if the file already exists
        if os.path.exists(h5_filename):
            print(f'Skipping computation for Chromosome {chromosome}, Resolution {resolution}: File already exists.')
            continue
        
        # Get the dark bins
        try:
            dark_bins = get_dark_bins(chromosome, resolution, bw)
        except ValueError as e:
            print(e)
            continue
        
        # Convert the indices to a NumPy array
        result_array = np.array(dark_bins)
        
        # Save the indices to an H5 file
        with h5py.File(h5_filename, 'w') as h5file:
            h5file.create_dataset('dark_bins_indices', data=result_array)
        
        # Print the number of removed bins
        print(f'Chromosome {chromosome}, Resolution {resolution}: {len(dark_bins)} bins removed.')

# Close the BigWig file
bw.close()

