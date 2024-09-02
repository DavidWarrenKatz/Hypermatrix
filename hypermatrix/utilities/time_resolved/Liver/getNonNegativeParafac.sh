#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=3-9:00:00
#SBATCH --ntasks=5
#SBATCH --mem=150000
#SBATCH --error=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetStructuredDataLog.err
#SBATCH --output=/home/dwk681/workspace/CRA004660/Liver/logs/%JgetStructuredDataLog.out


module load matlab/r2022b

# Set the path for data storage
data_path="/home/dwk681/workspace/CRA004660/Liver/"

# Start of MATLAB code
matlab -nodisplay -nosplash << MATLAB_SCRIPT
path = '$data_path';
addpath(genpath('/projects/b1198/epifluidlab/david'));


% Step 1: Load the HDF5 file
file_path = '/home/dwk681/workspace/CRA004660/Brain/tensor_normalized_filtered_106_500_32285.h5';
info = h5info(file_path);

% Step 2: Load the tensor
tensor = h5read(file_path, '/tensor');

% Step 3: Load the labels
labels = h5read(file_path, '/labels');

% Step 4: Convert labels back to a cell array of strings
labels = cellstr(labels);

        r = 3;
        sol_cell_downsampled = cell(1,r);
        sol_factors_downsampled = cell(1,r);
        ranks = 1:r;

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
            model.factorizations.myfac.data = T;
            model.factorizations.myfac.cpd  = {'U', 'V', 'S'};
            options.Display = 100;
            %options.TolFun=1e-50;
            %options.TolX=1e-50;
            options.MaxIter=300;
            sol = sdf_nls(model,options);
            sol_factors{1,i} = {sol.factors.U, sol.factors.V, sol.factors.S};
        end
        % Save the processed data for each chromosome
        output_file = sprintf('%sWorkspaces/structuredData_first3_300iterations.mat', path);
        save(output_file);
MATLAB_SCRIPT


