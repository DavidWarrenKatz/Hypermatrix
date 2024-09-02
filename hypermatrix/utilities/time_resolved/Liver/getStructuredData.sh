#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --time=3-9:00:00
#SBATCH --ntasks=5
#SBATCH --mem=150000
#SBATCH --error=/home/dwk681/workspace/CRA004660_new/Liver/logs/%JgetStructuredDataLog.err
#SBATCH --output=/home/dwk681/workspace/CRA004660_new/Liver/logs/%JgetStructuredDataLog.out


module load matlab/r2022b

# Set the path for data storage
data_path="/home/dwk681/workspace/CRA004660_new/Liver/"

# Start of MATLAB code
matlab -nodisplay -nosplash << MATLAB_SCRIPT
path = '$data_path';
addpath(genpath('/projects/b1198/epifluidlab/david'));

        % Load the correlation matrix
        file_path = sprintf('%sall_four_genes_filtered_count_matrices.mat', path);
        % Read the dataset from the HDF5 file
        A1 = h5read(file_path, '/Iso_Y');
        A2 = h5read(file_path, '/Iso_O');
        A3 = h5read(file_path, '/Het_Y');
        A4 = h5read(file_path, '/Het_O');

        r = 3;
        sol_cell_downsampled = cell(1,r);
        sol_factors_downsampled = cell(1,r);
        ranks = 1:r;

        n = size(A1,1);
        m = size(A1,2);
        size_tensor = [n m 4];
        R = 2;
        U = cpd_rnd(size_tensor,R);
        T = cpdgen(U);

        T(:,:,1) = A1;
        T(:,:,2) = A2;
        T(:,:,3) = A3;
        T(:,:,4) = A4;

        for i = 1:r
            model = struct;
            model.variables.u = randn(n,i);
            model.variables.v = randn(m,i);
            model.variables.s = randn(4,i);
            model.factors.U = {'u',@struct_nonneg};
            model.factors.V = {'v',@struct_nonneg};
            model.factors.S = {'s',@struct_nonneg};
            model.factorizations.myfac.data = T;
            model.factorizations.myfac.cpd  = {'U', 'V', 'S'};
            options.Display = 100;
            %options.TolFun=1e-50;
            %options.TolX=1e-50;
            options.MaxIter=500;
            sol = sdf_nls(model,options);
            sol_factors{1,i} = {sol.factors.U, sol.factors.V, sol.factors.S};
        end
        % Save the processed data for each chromosome
        output_file = sprintf('%sWorkspaces/minRows_genefiltered_structuredData_all4_first3_500iterations_highmem.mat', path);
        save(output_file);
MATLAB_SCRIPT


