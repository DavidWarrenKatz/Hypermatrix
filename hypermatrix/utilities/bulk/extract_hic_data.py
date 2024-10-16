import os
import sys
import numpy as np
import hicstraw
import h5py

def extract_hic_data(path, hic_url, resolutions, chromosomes, data_types):
    hic = hicstraw.HiCFile(hic_url)
    normalization = 'KR'
    resolutions = [int(r.strip().strip("'")) for r in resolutions.split(",")]
    chromosomes = [c.strip().strip("'") for c in chromosomes.split(",")]
    data_types = [d.strip().strip("'") for d in data_types.split(",")]

    for resolution in resolutions:
        for chromosome in chromosomes:
            startingbp = 0
            endingbp = hic.getChromosomes()[int(chromosome)].length
            for data in data_types:
                matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, normalization, "BP", resolution)
                matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                
                # Save short score to text file
                short_score_path = f'{path}hicFiles/short_score_textform/shortScore_res{resolution}_ch{chromosome}_{data}_{normalization}.txt'
                with open(short_score_path, 'wt') as fid_org:
                    for i in range(matrix.shape[0]):
                        for j in range(i, matrix.shape[1]):
                            if matrix[i, j] != 0:
                                fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(
                                    chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))
                print(f"[DEBUG] Saved short score file to: {short_score_path}")

                # Save matrix directly as .h5 file
                h5_filename = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data}_{normalization}.h5'
                with h5py.File(h5_filename, 'w') as hf:
                    hf.create_dataset('matrix', data=matrix)
                print(f"[DEBUG] Saved .h5 file to: {h5_filename}")

    # Generate combined short score files for all chromosomes
    for data in data_types:
        combined_short_score_path = f'{path}hicFiles/short_score_textform/shortScore_all_{data}_KR.txt'
        with open(combined_short_score_path, 'wt') as fid_org:
            for resolution in resolutions:
                for chromosome in chromosomes:
                    startingbp = 0
                    endingbp = hic.getChromosomes()[int(chromosome)].length
                    matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, "KR", "BP", resolution)
                    matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                    for i in range(matrix.shape[0]):
                        for j in range(i, matrix.shape[1]):
                            if matrix[i, j] != 0:
                                fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(
                                    chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))
        print(f"[DEBUG] Saved combined short score file to: {combined_short_score_path}")

if __name__ == "__main__":
    data_path = sys.argv[1]
    hic_url = sys.argv[2]
    resolutions_list = sys.argv[3]
    chromosomes_list = sys.argv[4]
    data_types_list = sys.argv[5]

    print(f"[DEBUG] Input data_path: {data_path}")
    print(f"[DEBUG] Input hic_url: {hic_url}")
    print(f"[DEBUG] Input resolutions_list: {resolutions_list}")
    print(f"[DEBUG] Input chromosomes_list: {chromosomes_list}")
    print(f"[DEBUG] Input data_types_list: {data_types_list}")

    extract_hic_data(data_path, hic_url, resolutions_list, chromosomes_list, data_types_list)
