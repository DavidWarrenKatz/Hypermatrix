import sys
import numpy as np
import hicstraw
import h5py

def extract_hic_data(path, hic_url, resolutions, chromosomes, data_types):
    print(f"[DEBUG] Starting extract_hic_data with hic_url: {hic_url}")
    try:
        hic = hicstraw.HiCFile(hic_url)
        print(f"[DEBUG] Successfully created HiCFile object")
    except Exception as e:
        print(f"[ERROR] Failed to create HiCFile object: {e}")
        return

    normalization = 'KR'
    try:
        resolutions = [int(r.strip().strip("'")) for r in resolutions.split(",")]
        chromosomes = [c.strip().strip("'") for c in chromosomes.split(",")]
        data_types = [d.strip().strip("'") for d in data_types.split(",")]
        print(f"[DEBUG] Parsed resolutions: {resolutions}")
        print(f"[DEBUG] Parsed chromosomes: {chromosomes}")
        print(f"[DEBUG] Parsed data_types: {data_types}")
    except Exception as e:
        print(f"[ERROR] Failed to parse inputs: {e}")
        return

    for resolution in resolutions:
        for chromosome in chromosomes:
            startingbp = 0
            try:
                endingbp = hic.getChromosomes()[int(chromosome)].length
                print(f"[DEBUG] Chromosome {chromosome} length: {endingbp}")
            except Exception as e:
                print(f"[ERROR] Failed to get chromosome length: {e}")
                continue

            for data in data_types:
                try:
                    matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, normalization, "BP", resolution)
                    matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                    print(f"[DEBUG] Retrieved matrix for chromosome {chromosome}, resolution {resolution}, data {data}")
                except Exception as e:
                    print(f"[ERROR] Failed to get matrix data: {e}")
                    continue

                try:
                    with open(f'{path}hicFiles/short_score_textform/shortScore_res{resolution}_ch{chromosome}_{data}_{normalization}.txt', 'wt') as fid_org:
                        for i in range(matrix.shape[0]):
                            for j in range(i, matrix.shape[1]):
                                if matrix[i, j] != 0:
                                    fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))
                    print(f"[DEBUG] Saved short form text for chromosome {chromosome}, resolution {resolution}, data {data}")
                except Exception as e:
                    print(f"[ERROR] Failed to save short form text: {e}")

                try:
                    new_filename = f'{path}Workspaces/individual/ch{chromosome}_res{resolution}_{data}_{normalization}.mat'
                    with h5py.File(new_filename, 'w') as hf:
                        hf.create_dataset('matrix', data=matrix)
                    print(f"[DEBUG] Saved .mat file for chromosome {chromosome}, resolution {resolution}, data {data}")
                except Exception as e:
                    print(f"[ERROR] Failed to save .mat file: {e}")

    for data in data_types:
        try:
            with open(f'{path}hicFiles/short_score_textform/shortScore_all_{data}_KR.txt', 'wt') as fid_org:
                for resolution in resolutions:
                    for chromosome in chromosomes:
                        startingbp = 0
                        endingbp = hic.getChromosomes()[int(chromosome)].length
                        matrix_object = hic.getMatrixZoomData(chromosome, chromosome, data, "KR", "BP", resolution)
                        matrix = matrix_object.getRecordsAsMatrix(startingbp, endingbp, startingbp, endingbp)
                        for i in range(matrix.shape[0]):
                            for j in range(i, matrix.shape[1]):
                                if matrix[i, j] != 0:
                                    fid_org.write('0 {} {} 0 0 {} {} 1 {:.5f}\n'.format(chromosome, startingbp + resolution * i, chromosome, startingbp + resolution * j, matrix[i, j]))
            print(f"[DEBUG] Saved combined short form text for data {data}")
        except Exception as e:
            print(f"[ERROR] Failed to save combined short form text: {e}")

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
