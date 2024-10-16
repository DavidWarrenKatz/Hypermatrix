addpath(genpath('../../bin/softwarefiles'));

% Load configuration file
config = load('../../utilities/single-cell/config.mat');
output_directory = config.output_directory;
resolutions = config.resolutions;
iterations = config.iterations;
filtered_list = config.filtered_list;

% Convert resolutions to cell array of strings if they are not already
if ischar(resolutions)
    resolutions = {resolutions};
elseif isstring(resolutions)
    resolutions = cellstr(resolutions);
end

% Initialize cell arrays to store resolutions and labels
resolution_values = cell(1, numel(resolutions));
resolution_labels = cell(1, numel(resolutions));

% Use regexp to find the integer before the colon and the label after the colon
for i = 1:numel(resolutions)
    tokens = regexp(resolutions{i}, '(\d+):(\w+)', 'tokens');
    if ~isempty(tokens)
        resolution_values{i} = str2double(tokens{1}{1});
        resolution_labels{i} = tokens{1}{2};
    end
end

% Display the extracted resolutions and labels
disp('Resolutions:');
disp(resolution_values);
disp('Labels:');
disp(resolution_labels);

% Define the prefix list file path
prefix_file_path = filtered_list;

% Read prefixes from the file
fileID = fopen(prefix_file_path, 'r');
if fileID == -1
    error('Could not open file %s for reading.', prefix_file_path);
end
prefixes = textscan(fileID, '%s');
fclose(fileID);
prefixes = prefixes{1};

% Define the list of chromosomes
chromosomes = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'};

for j = 1:numel(chromosomes)
    chromosome = chromosomes{j};
    for k = 1:numel(prefixes)
        prefix = prefixes{k};
        
        % Iterate over all resolutions
        for l = 1:numel(resolution_values)
            resolution = resolution_values{l};
            label = resolution_labels{l};
            
            % Construct file paths
            tensor_file_path = sprintf('%s/%s_combined_cumulant/%s/%s_%s_combined_cumulant.h5', output_directory, label, chromosome, prefix, chromosome);
            output_file_U = sprintf('%s/tensor_%s_AB_factors_cumulant/%s/%s_weights.h5', output_directory, label, chromosome, prefix);
            output_file_V = sprintf('%s/tensor_%s_AB_factors_cumulant/%s/%s_compartments.h5', output_directory, label, chromosome, prefix);

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
            
            % Check if the tensor file exists before trying to load it
            if ~exist(tensor_file_path, 'file')
                fprintf('Tensor file %s does not exist. Skipping.\n', tensor_file_path);
                continue;
            end

            try
                % Load the combined cumulant tensor from the HDF5 file
                combined_tensor = h5read(tensor_file_path, '/combined_cumulant_tensor');
                tensor_size = size(combined_tensor);
                fprintf('Processing combined tensor of size %s for %s.\n', mat2str(tensor_size), tensor_file_path);

                % Perform the low-rank tensor factorization
                model = struct;
                model.variables.u = randn(size(combined_tensor, 1), 2);  % Assuming rank 2
                model.variables.v = randn(size(combined_tensor, 3), 2);
                model.factors.U = {'u', @struct_nonneg};
                model.factors.V = {'v', @struct_nonneg};
                model.factorizations.myfac.data = combined_tensor;
                model.factorizations.myfac.cpd = {'U', 'V', 'V', 'V'};
                options.Display = 100;
                options.MaxIter = iterations;
                options.TolFun = 1e-12;  % Set the tolerance here
                options.TolX = 1e-12;
                sol = sdf_nls(model, options);

                % Extract the required fields from the nested structure
                U = sol.factors.U;
                V = sol.factors.V;
                
                % Create and write to the HDF5 files
                output_dataset_U = '/weights';
                output_dataset_V = '/compartment_factors';
                
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
end

