import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.io import savemat
import math
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

window = "175kb"
path="/projects/b1198/epifluidlab/david/GSE63525/GM12878/"
resolutions = [25000]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 20', '21', 22']
threshold = .8

# Define a function to write valley crossings to a BED file
def save_bed_file(valley_crossings, filename, ch):
    with open(filename, 'w') as bed_file:
        for index in valley_crossings.keys():
         # Create BED format entry (chromosome, start, end)
         bed_entry = f'{ch}\t{index}\t{index+1}\t.\t{valley_crossings[index]}\t+\n'
         bed_file.write(bed_entry)

for resolution in resolutions:
    for ch in chromosomes:
        mat_file = "/projects/b1198/epifluidlab/david/GSE63525/GM12878/Workspaces/individual/ch{ch}_res{resolution}_structedData_10Downsampled_oe_log2_KR_first4_200iterations.mat"
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
        max_hypermatrix_smoothing_factor = None
        max_insulation_smoothing_factor = None

        removed_rows_path = path + f'Workspaces/individual/ch{ch}_res{resolution}_oe_KR_removedRows.mat'    
        data = scipy.io.loadmat(removed_rows_path)
        removed_indices = data['removed_rows_indices']
        # Convert the indices to a NumPy array
        result_array = np.array(removed_indices)
        removed_rowIndices_by_ch[ch] = result_array
        N_ch = factor_one.size
        filtered_indices_ch = [i for i in range(N_ch-1) if i not in removed_rowIndices_by_ch[ch]]

        for idx, vector in enumerate([vector1, vector2, vector3, vector4]):
            # Extract the insulation vector from the second set
            insulation_vector = insulation_by_window_ch[window][ch]
            # Define the threshold below which values will be set
            threshold = -5  # You can adjust the threshold as needed
            # Set values below the threshold to be equal to the threshold
            insulation_vector[insulation_vector < threshold] = threshold
            # Define the threshold below which values will be set
            positive_threshold = 40  # You can adjust the threshold as needed
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
                max_hypermatrix_smoothing_factor = hypermatrix_smoothing_factor
                max_insulation_smoothing_factor = insulation_vector_smoothing_factor

        # Define your curve as a numpy array
        curve = max_corr_vector_pair[0] 
        
        local_minima_indices = argrelextrema(curve, np.less)[0]

        # Calculate the delta vector
        delta_vector = np.diff(curve)

        # Find zero-crossings in the delta vector
        zero_crossings = np.where(np.diff(np.sign(delta_vector)))[0]

        # Filter out zero-crossings occurring at peaks
        valley_crossings = {}
        for crossing in zero_crossings:
            left_max = max(delta_vector[:crossing], default=0)
            right_min = min(delta_vector[crossing+1:], default=0)
            boundary_strength = left_max - right_min
            if boundary_strength >= threshold:
                valley_crossings[crossing] = boundary_strength

        # Create a Venn diagram to visualize the common local minima
        venn2([set(valley_crossings.keys()), set(local_minima_indices)], ('Valley Crossings', 'Local Minima'))
         
       # Save the Venn diagram figure to a file
        venn_diagram_filename = 'venn_diagram.png'
        plt.savefig(venn_diagram_filename)
        plt.show()

        plt.savefig(path + f'Images/res_{resolution}_ch{ch}_min_venn.png')
        # Close the figure 
        plt.close(fig)

         # Save valley crossings to a BED file
         bed_filename = path + 'boundaries/filtered_valleys_res{resolution}_ch{ch}_oe_log2_KR_threshold_{threshold}.bed'
         save_bed_file(valley_crossings, bed_filename)
         print(f'Filtered valleys saved to {bed_filename}')

