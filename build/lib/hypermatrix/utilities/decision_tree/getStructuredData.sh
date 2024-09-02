#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomicslong
#SBATCH --time=6-15:00:00
#SBATCH --ntasks=5
#SBATCH --mem=150000
#SBATCH --error=/home/dwk681/workspace/GSE141252_Stoeger_et_al._2022/logs/%JgetStructuredDataLog.err
#SBATCH --output=/home/dwk681/workspace/GSE141252_Stoeger_et_al._2022/logs/%JgetStructuredDataLog.out


module load matlab/r2022b

# Set the path for data storage
data_path="/home/dwk681/workspace/GSE141252_Stoeger_et_al._2022"

# Start of MATLAB code
matlab -nodisplay -nosplash << MATLAB_SCRIPT
% MATLAB code here
path = '$data_path';
addpath(genpath('/projects/b1198/epifluidlab/david'));
        
        filename = '/home/dwk681/workspace/GSE141252_Stoeger_et_al._2022/tensor_data/RNA_tensor.h5'
        T = h5read(filename, '/zscored_tensor');

        r = 100;
        sol_cell_downsampled = cell(1,r);
        sol_factors_downsampled = cell(1,r);
        ranks = 1:r;

        n = size(T, 1);
        m = size(T, 2);
        p = size(T, 3);

     % Factorization loop
       i = 100
         model = struct;
         model.variables.u = randn(n, i);
         model.variables.v = randn(m, i);
         model.variables.p = randn(p, i); 
         model.factors.U = {'u', @struct_nonneg};
         model.factors.V = {'v', @struct_nonneg};
         model.factors.P = {'p', @struct_nonneg}; 
         model.factorizations.myfac.data = T;
         model.factorizations.myfac.cpd = {'U', 'V', 'P'};
         options.Display = 100;
         options.MaxIter = 300;
         sol = sdf_nls(model, options);
         sol_factors{1, i} = {sol.factors.U, sol.factors.V, sol.factors.P};
    

% Save the processed data
output_file = sprintf('%s/tensor_data/zscored_tensor_rank20_300iterations.mat', path);
save(output_file, 'sol_factors');

MATLAB_SCRIPT



