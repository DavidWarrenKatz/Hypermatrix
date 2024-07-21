import sys
import numpy as np
import h5py

def read_matrix_from_text_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    matrix = []
    for line in lines:
        row = [float(val) if val != 'nan' else np.nan for val in line.strip().split()]
        matrix.append(row)

    matrix = np.array(matrix)
    return matrix

def process_pearsons(path, resolutions, chromosomes):
    resolutions = [int(r.strip()) for r in resolutions.split(",")]
    chromosomes = [c.strip().strip("'") for c in chromosomes.split(",")]

    for chromosome in chromosomes:
        for resolution in resolutions:
            filename = f"{path}pearsons/individual/res{resolution}_ch{chromosome}_observed_KR_pearsons.txt"
            matrix = read_matrix_from_text_file(filename)

            translated_matrix = matrix + 1

            new_filename = f"{path}Workspaces/individual/translated_ch{chromosome}_res{resolution}_pearsons_KR.mat"
            new_filename_npy = f"{path}Workspaces/individual/translated_ch{chromosome}_res{resolution}_pearsons_KR.npy"

            with h5py.File(new_filename, 'w') as hf:
                hf.create_dataset('translated_matrix', data=translated_matrix)

            np.save(new_filename_npy, translated_matrix)

if __name__ == "__main__":
    data_path = sys.argv[1]
    resolutions_list = sys.argv[2]
    chromosomes_list = sys.argv[3]
    process_pearsons(data_path, resolutions_list, chromosomes_list)
