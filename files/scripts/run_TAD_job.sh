#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=02:50:00
#SBATCH --ntasks=5
#SBATCH --mem=100000
#SBATCH --error=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JFindTADS.err
#SBATCH --output=/projects/b1198/epifluidlab/david/GSE63525/GM12878/logs/%JFINDTADS.out


python<<EOF
#This code makes the images for TAD calling comparisons

import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import scipy
from scipy.io import savemat
import math
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import textwrap
from scipy.stats import pearsonr, spearmanr
import math
from scipy.stats import t
from scipy.signal import convolve
import pyBigWig


window = "175kb"
path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"
grant_path="/home/dwk681/workspace/grant/"
resolutions = [25000]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

removed_rowIndices_by_ch = {}
insulation_by_window_ch = {}

# Define a function to write valley crossings to a BED file
def save_bed_file(valley_crossings, filename, ch, resolution, scanning_window):
    with open(filename, 'w') as bed_file:
        for index in valley_crossings.keys():
            # Create BED format entry (chromosome, start, end)
            start_position = int(index) * resolution 
            end_position = (int(index) + 1) * resolution 
            bed_entry = f'{ch}	{start_position}	{end_position}	.	{valley_crossings[index]}	+
'
            bed_file.write(bed_entry)

def find_local_extrema(vector):
    extrema = {}  # Initialize an empty dictionary to store extrema

    for i in range(1, len(vector) - 1):
        if vector[i] < vector[i - 1] and vector[i] < vector[i + 1]:
            extrema[i] = "min"  # Local minimum
        elif vector[i] > vector[i - 1] and vector[i] > vector[i + 1]:
            extrema[i] = "max"  # Local maximum

    return extrema

def has_previous_local_max(index, extrema):
    previous_max_found = False
    index_keys = list(extrema.keys())
    index_keys_before = [i for i in index_keys if i < index]
    for extrema_index in index_keys_before:
        if extrema[extrema_index] == "max":
            previous_max_found = True
            break
    return previous_max_found

def has_local_min_to_the_right(index, extrema):
    local_min_found = False
    index_keys = list(extrema.keys())
    for extrema_index in index_keys:
        if extrema_index <= index:
            continue
        if extrema[extrema_index] == "min":
            local_min_found = True
            break
    return local_min_found

def find_local_min_to_the_right(index, extrema):
    index_keys = list(extrema.keys())
    for extrema_index in index_keys:
        if extrema_index <= index:
            continue
        if extrema[extrema_index] == "min":
            return extrema_index
    return None  # Return None if no local minimum to the right is found

def find_local_max_to_the_left(index, extrema):
    index_keys = list(extrema.keys())
    local_max_index = None

    for extrema_index in reversed(index_keys):
        if extrema_index >= index:
            continue
        if extrema[extrema_index] == "max":
            local_max_index = extrema_index
            break

    return local_max_index

def filter_local_mins(extrema, vector, threshold):
    filtered_local_mins = {}
    
    for index, extrema_type in extrema.items():
        if extrema_type == "min":
            left_max = find_local_max_to_the_left(index, extrema)
            right_min = find_local_min_to_the_right(index, extrema)
            if left_max==None or right_min==None:
                filtered_local_mins[index] = 0
            else:
                if left_max - right_min > threshold:
                    filtered_local_mins[index] = left_max - right_min
   
    return filtered_local_mins

def save_curve_as_bigwig(curve, ch, resolution, output_path):
    # Create a new BigWig file
    bw = pyBigWig.open(output_path, "w")

    # Prepare data in the required format for writing to the BigWig file
    chrom = ch
    start = 0  # Start position for the chromosome
    step = resolution  # The resolution should be used as the step size
    span = resolution  # The span is also set to the resolution
    values = curve  # The values for the curve

    # Write data to the BigWig file
    bw.addEntries(chrom, start, ends=None, values=values, span=span, step=step)

    # Close the BigWig file
    bw.close()


for resolution in resolutions:
    for ch in chromosomes:
        mat_file = f"/projects/b1198/epifluidlab/david/GSE63525/GM12878/Workspaces/individual/ch{ch}_res{resolution}_structedData_10Downsampled_oe_log2_KR_first4_200iterations.mat"
        data = scipy.io.loadmat(mat_file)
        factor_one = data['sol_factors_downsampled'][0, 3][0, 0][:, 0]
        factor_two = data['sol_factors_downsampled'][0, 3][0, 0][:, 1]
        factor_three = data['sol_factors_downsampled'][0, 3][0, 0][:, 2]
        factor_four = data['sol_factors_downsampled'][0, 3][0, 0][:, 3]

        # Initialize variables to store the most correlated vectors and their correlation values
        max_corr_vector_pair = None
        max_corr_value = -1
        max_corr_window = None
        max_corr_hypermatrix = None

        loaded_data = np.load(path + f'Images/dataForInsulation_res{resolution}_oe_KR_log2.npy', allow_pickle=True).item()
        insulation_by_window_ch = loaded_data['insulation_by_window_ch']

        removed_rows_path = path + f'Workspaces/individual/ch{ch}_res{resolution}_oe_KR_removedRows.mat'
        data = scipy.io.loadmat(removed_rows_path)
        removed_indices = data['removed_rows_indices']
        # Convert the indices to a NumPy array
        result_array = np.array(removed_indices)
        removed_rowIndices_by_ch[ch] = result_array
        N_ch = factor_one.size
        filtered_indices_ch = [i for i in range(N_ch-1) if i not in removed_rowIndices_by_ch[ch]]

        for idx, vector in enumerate([factor_one, factor_two, factor_three, factor_four]):
            # Extract the insulation vector from the second set
            insulation_vector = insulation_by_window_ch[window][ch]
            # Define the threshold below which values will be set
            threshold = -5  # You can adjust the threshold as needed
            # Set values below the threshold to be equal to the threshold
            insulation_vector[insulation_vector < threshold] = threshold
            # Define the threshold below which values will be set
            positive_threshold = 30  # You can adjust the threshold as needed
            # Set values below the threshold to be equal to the threshold
            vector[vector > positive_threshold] = threshold
            # Calculate Pearson correlation coefficient and p-value for smoothed vectors
            corr_value, p_value = pearsonr(vector[filtered_indices_ch], insulation_vector[filtered_indices_ch])

            # Update best correlation results if a better correlation is found
            if corr_value > max_corr_value:
                max_corr_value = corr_value
                max_corr_window = window
                max_corr_hypermatrix = idx + 1
                max_corr_vector_pair = (vector, insulation_vector)

        # Define your curve as a numpy array
        curve = max_corr_vector_pair[0]
        
        output_path = grant_path + f'Worksapces/bigWigs/ch{ch}_res{resolution}_Hypermatrix_Tad_10Downsampled_oe_log2_KR_first4_200iterations.bigwig'
        save_curve_as_bigwig(curve, ch, resolution, output_path)
        
        local_minima_indices = argrelextrema(curve, np.less)[0]

        # Set the threshold for the difference between local max and local min
        threshold = .2  # Adjust this threshold as needed
        scanning_window = 5;

        filtered_minima = {}

        # Iterate through the local minima and apply the threshold condition
        for min_index in local_minima_indices:  
            left_max = max(curve[min_index - scanning_window: min_index]) if min_index > scanning_window else curve[min_index]
            right_min = min(curve[min_index:min_index+scanning_window]) if min_index < len(curve) - 1 else curve[min_index]

            if left_max - right_min >= threshold:
                filtered_minima[min_index] = left_max - right_min

        threshold_new = .001        
        extrema = find_local_extrema(curve)
        filtered_extrema = filter_local_mins(extrema, curve, threshold_new)    
            
        # Create a Venn diagram to visualize the common local minima
        venn2([set(filtered_minima), set(filtered_extrema)], ('Filtered Minima', 'Filtered Extrema'))

        # Save the Venn diagram figure to a file
        venn_diagram_filename = 'venn_diagram.png'
        plt.savefig(venn_diagram_filename)
        plt.show()

        plt.savefig(grant_path + f'Images/res_{resolution}_ch{ch}_min_venn.png')

        # Save valley crossings to a BED file
        bed_filename = grant_path + f'Workspaces/boundaries/filtered_valleys_res{resolution}_ch{ch}_oe_log2_KR_threshold_{threshold}_scanning_{scanning_window}.bed'
        save_bed_file(filtered_minima, bed_filename, ch, resolution, scanning_window)
        print(f'Filtered valleys saved to {bed_filename}')
EOF
