#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=10:00:00
#SBATCH --ntasks=5
#SBATCH --mem=150000
#SBATCH --error=/projects/b1198/epifluidlab/david/grant/logs/%JgetStructuredData.err
#SBATCH --output=/projects/b1198/epifluidlab/david/grant/logs/%JgetStructuredData.out

# Load the necessary modules for conda
source /etc/profile
source ~/.bashrc

# Activate the conda environment
conda activate multiomics6

module load matlab/r2022b

# Set the path for data storage
data_path="/home/dwk681/workspace/grant/"

# Start of MATLAB code
matlab<<MATLAB_SCRIPT
% MATLAB code here
path = '$data_path';
addpath(genpath('/home/dwk681/workspace/softwareFiles/tensorlab'));

resolutions = [1000000];
chromosomes = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'};
for resolution = resolutions
    for j = 1:numel(chromosomes)
        chromosome = chromosomes{j};
        file_path = sprintf('%sWorkspaces/mCpG_chr%s_res1000000_small_10000_ks_p_value_correlation.h5', path, chromosome);
        dataset_name = '/matrix';  % Name of the dataset within the HDF5 file;
        % Read the dataset from the HDF5 file
        matrix = h5read(file_path, dataset_name);
        
        A1 = matrix;
        A2 = downsample_matrix(matrix, 50);

        r = 3;
        sol_cell = cell(1,r);
        sol_factors = cell(1,r);
        ranks = 1:r;

        n = size(A1,1);
        size_tensor = [n n 2];
        R = 2;
        U = cpd_rnd(size_tensor,R);
        T = cpdgen(U);

        T(:,:,1) = A1;
        T(:,:,2) = A2;

        for i = 1:r
            model = struct;
            model.variables.u = randn(n,i);
            model.variables.s = randn(2,i);
            model.factors.U = {'u',@struct_nonneg};
            model.factors.S = {'s',@struct_nonneg};
            model.factorizations.myfac.data = T;
            model.factorizations.myfac.cpd  = {'U', 'U', 'S'};
            options.Display = 100;
            %options.TolFun=1e-50;
            %options.TolX=1e-50;
            options.MaxIter=400;
            sol = sdf_nls(model,options);
            sol_factors{1,i} = {sol.factors.U, sol.factors.S};
        end
        % Save the processed data for each chromosome
        output_file = sprintf('%sWorkspaces/individual/ch%s_res%d_structedData_2Downsampled_mCpG_small_10000_first3_400iterations.mat', path, chromosome, resolution);
        save(output_file);
     end
end

MATLAB_SCRIPT
