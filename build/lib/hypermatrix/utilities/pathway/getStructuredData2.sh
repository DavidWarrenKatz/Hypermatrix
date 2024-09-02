#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=5:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100000
#SBATCH --error=/home/dwk681/workspace/GSE143519/logs/%JgetStructuredDataLog.err
#SBATCH --output=/home/dwk681/workspace/GSE143519/logs/%JgetStructuredDataLog.out


module load matlab/r2022b

# Set the path for data storage
data_path="/home/dwk681/workspace/GSE143519/"

# Start of MATLAB code
matlab -nodisplay -nosplash << MATLAB_SCRIPT
path = '$data_path';
addpath(genpath('/projects/b1198/epifluidlab/david'));

        
        filename = '/home/dwk681/workspace/GSE143519/data/RNA_massSpec_data_823_2_2_2_tensor.h5'
        T = h5read(filename, '/tensor');

        r = 10;
        sol_cell_downsampled = cell(1,r);
        sol_factors_downsampled = cell(1,r);
        ranks = 1:r;


        n = size(T, 1);
        m = size(T, 2);
        p = size(T, 3);
        q = size(T, 4);

        for i = 1:r
            model = struct;
            model.variables.u = randn(n,i);
            model.variables.v = randn(m,i);
            model.variables.s = randn(p,i);
            model.variables.w = randn(q,i);

            model.factors.U = {'u',@struct_nonneg};
            model.factors.V = {'v',@struct_nonneg};
            model.factors.S = {'s',@struct_nonneg};
            model.factors.W = {'w',@struct_nonneg};

            model.factorizations.myfac.data = T;
            model.factorizations.myfac.cpd  = {'U', 'V', 'S', 'W'};
            options.Display = 100;
            %options.TolFun=1e-50;
            %options.TolX=1e-50;
            options.MaxIter=800;
            sol = sdf_nls(model,options);
            sol_factors{1,i} = {sol.factors.U, sol.factors.V, sol.factors.S, sol.factors.W};
        end
        % Save the processed data for each chromosome
        output_file = sprintf('%sWorkspaces/translated_structuredData_massSpec_RNA_first10_800iterations_zscores.mat', path);
        save(output_file);
MATLAB_SCRIPT

