% getStructuredData.m
addpath(genpath('../projects/softwarefiles/tensorlab'));

% Load configuration file
config = load('config.mat');
data_path = config.data_path;
resolutions = config.resolutions;
chromosomes = cellstr(config.chromosomes); % Convert to cell array of character vectors if needed
iterations = config.iterations;

for resolution = resolutions
    for j = 1:numel(chromosomes)
        chromosome = chromosomes{j};
        file_path = sprintf('%sWorkspaces/individual/ch%s_res%d_oe_KR_cumulant.h5', data_path, chromosome, resolution);
        output_file = sprintf('%sWorkspaces/individual/ch%s_res%d_structuredData_3rdCumulant_rank2_%diterations.h5', data_path, chromosome, resolution, iterations);
        
        % Check if the output file already exists
        if exist(output_file, 'file') == 2
            fprintf('File %s already exists. Skipping computation.\n', output_file);
            continue;
        end
        
        dataset_name = '/degree_3_cumulant';
        cumulant = h5read(file_path, dataset_name);
        model = struct;
        model.variables.u = randn(size(cumulant, 1), 2);
        model.factors.U = {'u', @struct_nonneg};
        model.factorizations.myfac.data = cumulant;
        model.factorizations.myfac.cpd = {'U', 'U', 'U'};
        options.Display = 100;
        options.MaxIter = iterations;
        sol = sdf_nls(model, options);
        
        % Extract the required field from the nested structure
        U = sol.factors.U;
        
        % Create and write to the HDF5 file
        output_dataset = '/U';
        h5create(output_file, output_dataset, size(U));
        h5write(output_file, output_dataset, U);
    end
end

