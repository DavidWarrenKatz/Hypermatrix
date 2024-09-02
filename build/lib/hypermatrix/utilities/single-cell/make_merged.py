import numpy as np
import os
import argparse
from sklearn.preprocessing import quantile_transform

# Parse command-line arguments for resolution and label
parser = argparse.ArgumentParser(description='Merge compartment calls for different resolutions.')
parser.add_argument('--resolution', type=str, help='Resolution value, e.g., "1000000" for 1Mb')
parser.add_argument('--label', type=str, help='Resolution label, e.g., "1Mb"')
parser.add_argument('--indir_base', type=str, help='Base input directory, e.g., "/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/new_processing3/"')

args = parser.parse_args()

# Set variables based on command-line arguments
res0 = args.label
mode = 'pad1_std1_rp0.5_sqrtvc'
indir = f'{args.indir_base}hicluster_{res0}_impute_dir/merged/'

comp = []
for c in range(1, 23):
    tmp_file = f'{indir}{mode}_chr{c}.cpgcomp.npy'
    if os.path.exists(tmp_file):
        tmp = np.load(tmp_file)
        binfilter = (np.std(tmp, axis=0) > 0)
        comptmp = np.ones(tmp.shape) / 2
        comptmp[:, binfilter] = quantile_transform(tmp[:, binfilter], output_distribution='uniform', n_quantiles=int(np.sum(binfilter)//10), axis=1)
        comp.append(comptmp)
        print(c)

comp = np.concatenate(comp, axis=1)
output_file = f'{indir}all_merged_{res0}_cpgcmp.{mode}.txt'
np.savetxt(output_file, comp, fmt='%s', delimiter='\t')
