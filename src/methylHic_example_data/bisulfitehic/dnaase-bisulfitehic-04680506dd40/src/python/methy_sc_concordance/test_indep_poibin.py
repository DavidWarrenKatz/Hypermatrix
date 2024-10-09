##for each cpg, get mean and variance across alleles
#expected:cpg_mean=cpg1_mean+cpg2_mean, cpg_sd=cpg1_sd+cpg2_sd, then get the distribution
#obs: for each allele, must have two value in the same allele, get methy_sum, and calculate the mean and sd across alleles
#obs/exp >> 1 then is not

#input: a file contain all of cpg's methy value, get matrix in numpy array
#get colMenas, colVars and build distribution
#build large number of pairs in matrix, then expected sd: cpg1_sd+cpg2_sd; obs: sd(cpg1_methy+cpg2_methy)
import gzip,sys
import numpy as np


def test_indep_poibin(argv):
    input = argv[0]
    output = argv[1]
    max_dist = 3000000
    data_mat=[]
    cpg_pos=[]
    with gzip.open(input,'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            line = line.replace('NA','nan')
            #print(line + "\n")
            tmp_array = line.split("\t")
            array_row = np.array(tmp_array[6:],dtype='float32')
            data_mat.append(array_row)
            cpg_pos.append(tmp_array[2])
    data_mat = np.array(data_mat,dtype='float32')
    cpg_pos = np.array(cpg_pos,dtype='int32')
    #allele_mean = np.mean(data_mat,axis=1)
    #allele_std = np.std(data_mat, axis=1)
    with gzip.open(output,'wb') as out:
        for i in range(data_mat.shape[0]):
            cpgi = data_mat[i, :]
            cpgi_std = np.nanstd(cpgi)
            cpgi_pos = cpg_pos[i]
            ratio = []
            for j in range(data_mat.shape[0]):
                cpgj_pos = cpg_pos[j]
                #print(cpgj_pos)
                #print(cpgi_pos)
                if abs(cpgj_pos - cpgi_pos) > max_dist:
                    ratio.append('NA')
                cpgj = data_mat[j, :]
                cpgj_std = np.nanstd(cpgj)
                obs_std = np.nanstd((cpgi + cpgj) / 2)
                if not np.isnan(obs_std):
                    exp_std = cpgi_std + cpgj_std
                    ratio.append(obs_std/exp_std)
                else:
                    ratio.append('NA')
            outline = '\t'.join(np.array(ratio,dtype='str')) + '\n';
            out.write(outline.encode('utf-8'))


if __name__ == "__main__":
    test_indep_poibin(sys.argv[1:])
