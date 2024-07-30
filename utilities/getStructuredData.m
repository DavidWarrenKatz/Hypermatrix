% getStructuredData.m

data_path = '../projects/GSE63525/GM12878/';
addpath(genpath('../projects/softwarefiles/tensorlab'));

resolutions = [1000000];
chromosomes = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'};

for resolution = resolutions
    for j = 1:numel(chromosomes)
        chromosome = chromosomes{j};
        file_path = sprintf('%sWorkspaces/individual/ch%s_res%d_oe_KR_cumulant_new.h5', data_path, chromosome, resolution);
        dataset_name = '/degree_2_cumulant';
        cumulant = h5read(file_path, dataset_name);
        disp(['Size of cumulant tensor for chromosome ', chromosome, ': ', num2str(size(cumulant))]);

        file_path = sprintf('%sWorkspaces/individual/ch%s_res%d_oe_KR_nn_decomp.h5', data_path, chromosome, resolution);
        W = h5read(file_path, '/W');
        first_column_vector = W(:, 1);
        disp('Size of first column vector of the rank-2 decomposition:');
        disp(size(first_column_vector));

        r = 3;
        sol_cell = cell(1,r);
        sol_factors = cell(1,r);
        ranks = 1:r;
        n = size(cumulant,1);
        size_tensor = [n n n];
        R = 2;
        U = cpd_rnd(size_tensor,R);
        T = cpdgen(U);

        for i = 1:r
            model = struct;
            model.variables.u = [W, W];
            model.factors.U = {'u',@struct_nonneg};
            model.factorizations.myfac.data = T;
            model.factorizations.myfac.cpd = {'U', 'U', 'U'};
            options.Display = 100;
            options.MaxIter = 400;
            sol = sdf_nls(model,options);
            sol_factors{1,i} = {sol.factors.U};
        end

        output_file = sprintf('%sWorkspaces/individual/ch%s_res%d_structedData_2ndCumulant_first3_400iterations.mat', data_path, chromosome, resolution);
        save(output_file, 'sol_factors');
    end
end

