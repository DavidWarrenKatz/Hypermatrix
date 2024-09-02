#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=6-9:00:00
#SBATCH --ntasks=5
#SBATCH --mem=190000
#SBATCH --error=/home/dwk681/workspace/CRA004660/logs/%JgetStructuredDataLog.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/logs/%JgetStructuredDataLog.out


module load matlab/r2022b

# Set the path for data storage
data_path="/home/dwk681/workspace/CRA004660/"

# Start of MATLAB code
matlab -nodisplay -nosplash << MATLAB_SCRIPT
path = '$data_path';
addpath(genpath('/projects/b1198/epifluidlab/david'));

        
        file_path = sprintf('%snormalized_tensor.h5', path);
        tensor = h5read(file_path, '/tensor');

        r = 10;
        sol_factors = cell(1,r);

        n = size(tensor,1);
        m = size(tensor,2);
        p = size(tensor,3);

        for i = 1:r
            model = struct;
            model.variables.u = randn(n,i);
            model.variables.v = randn(m,i);
            model.variables.s = randn(p,i);
            model.factors.U = {'u',@struct_nonneg};
            model.factors.V = {'v',@struct_nonneg};
            model.factors.S = {'s',@struct_nonneg};
            model.factorizations.myfac.data = tensor;
            model.factorizations.myfac.cpd  = {'U', 'V', 'S'};
            options.Display = 100;
            %options.TolFun=1e-50;
            %options.TolX=1e-50;
            options.MaxIter=400;
            sol = sdf_nls(model,options);
            sol_factors{i} = {sol.factors.U, sol.factors.V, sol.factors.S};
        end
        
        output_file = sprintf('%sWorkspaces/structuredData_broken_first10_400iterations_highmem.mat', path);
        save(output_file);
MATLAB_SCRIPT


