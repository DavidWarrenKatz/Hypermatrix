import numpy as np
import pandas as pd
from config_and_print import data_path, resolutions, chromosomes

# Passed resolutions string from the shell script
resolutions_string = resolutions[0]
print(resolutions_string)

# Parse the resolutions string to a dictionary
resolutions = dict(res.split(":") for res in resolutions_string.split(","))

# Load chromosome sizes from a file specified by the chrom_file variable
chromsize = pd.read_csv('../../bin/softwarefiles/hg19.autosome.chrom.sizes', sep='\t', header=None, index_col=0).to_dict()[1]

# Loop over each resolution and chromosome to create BED segments
for res_str, label in resolutions.items():
    res = int(res_str)  # Convert resolution string to integer
    for c in chromsize:
        # Calculate the number of segments based on chromosome size and resolution
        ngene = int(chromsize[c] // res) + 1
        # Generate BED format data for each segment
        bed = [[c, i*res, (i+1)*res] for i in range(ngene)]
        # Ensure the last segment extends to the end of the chromosome
        bed[-1][-1] = chromsize[c]
        # Save the BED data to a file, one for each chromosome
        outdir = f'{output_directory}/hicluster_{label}_impute_dir/bins'
        # Make sure the output directory exists
        import os
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        np.savetxt(f'{outdir}/{c}.bed', bed, fmt='%s', delimiter='\t')
        print(f'Resolution: {res_str}, Chromosome: {c}')
