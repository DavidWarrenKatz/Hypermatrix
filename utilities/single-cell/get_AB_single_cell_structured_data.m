addpath(genpath('/home/dwk681/workspace/hypermatrix_test/hypermatrix/bin/softwarefiles'));

% Load configuration file
config = load('/home/dwk681/workspace/hypermatrix_test/hypermatrix/utilities/single-cell/config.mat');
data_path = config.output_directory;
resolutions = config.resolutions;
iterations = config.iterations;
filtered_list = config.filtered_list;

% Convert resolutions to cell array of strings if they are not already
if ischar(resolutions)
    resolutions = {resolutions};
elseif isstring(resolutions)
    resolutions = cellstr(resolutions);
end

% Use regexp to find the integer before the colon
tokens = regexp(resolutions, '(\d+):', 'tokens');

% Convert the result from cell array to integer
resolution = str2double(tokens{1}{1});

% Display the result
disp(resolution);

% Define the prefix list file path
prefix_file_path = '/home/dwk681/workspace/hypermatrix_test/hypermatrix/projects/single_cell_files/filtered_bam_list.txt';

% Read prefixes from the file
fileID = fopen(prefix_file_path, 'r');
prefixes = textscan(fileID, '%s');
fclose(fileID);
prefixes = prefixes{1};

% Define the list of chromosomes (you might need to adjust this part based on your actual chromosome list)
chromosomes = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'};

for j = 1:numel(chromosomes)
    chromosome = chromosomes{j};
    for k = 1:numel(prefixes)
        prefix = prefixes{k};
        
        % Construct file paths
        tensor_file_path = sprintf('/home/dwk681/workspace/hypermatrix_test/hypermatrix/projects/single_cell_files/hic_methy_1Mb_tensor_singlecell/%s/%s_%s.h5', chromosome, prefix, chromosome);
        output_file_U = sprintf('%s/tensor_AB_calls/%s/%s_res%s_structuredData_rank2_%diterations_U.h5', data_path, chromosome, resolution, iterations);
        output_file_V = sprintf('%s/tensor_AB_calls/%s/%s_res%s_structuredData_rank2_%diterations_V.h5', data_path, chromosome, resolution, iterations);

        % Ensure output directories exist
        [output_dir_U, ~, ~] = fileparts(output_file_U);
        [output_dir_V, ~, ~] = fileparts(output_file_V);
        if ~exist(output_dir_U, 'dir')
            mkdir(output_dir_U);
        end
        if ~exist(output_dir_V, 'dir')
            mkdir(output_dir_V);
        end

        % Check if the output files already exist
        if exist(output_file_U, 'file') == 2 && exist(output_file_V, 'file') == 2
            fprintf('Files %s and %s already exist. Skipping computation.\n', output_file_U, output_file_V);
            continue;
        end
        
        try
            % Load the tensor from the HDF5 file
            tensor = h5read(tensor_file_path, '/Tensor');
            tensor_size = size(tensor);
            fprintf('Processing tensor of size %s for %s.\n', mat2str(tensor_size), tensor_file_path);

            % Perform the tensor computation
            model = struct;
            model.variables.u = randn(size(tensor, 1), 2);
            model.variables.v = randn(size(tensor, 3), 2);
            model.factors.U = {'u', @struct_nonneg};
            model.factors.V = {'v', @struct_nonneg};
            model.factorizations.myfac.data = tensor;
            model.factorizations.myfac.cpd = {'U', 'V', 'V'};
            options.Display = 100;
            options.MaxIter = iterations;
            sol = sdf_nls(model, options);

            % Extract the required fields from the nested structure
            U = sol.factors.U;
            V = sol.factors.V;
            
            % Create and write to the HDF5 files
            output_dataset_U = '/U';
            output_dataset_V = '/V';
            
            h5create(output_file_U, output_dataset_U, size(U));
            h5write(output_file_U, output_dataset_U, U);
            
            h5create(output_file_V, output_dataset_V, size(V));
            h5write(output_file_V, output_dataset_V, V);
            
            fprintf('Computed and saved results to %s and %s.\n', output_file_U, output_file_V);
        catch ME
            fprintf('Error processing tensor for %s: %s\n', tensor_file_path, ME.message);
        end
    end
end
