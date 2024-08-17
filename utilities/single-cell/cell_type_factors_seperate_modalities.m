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

% Define the list of chromosomes (you might need to adjust this part based on your actual chromosome list)
chromosomes = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'};


for j = 1:numel(chromosomes)
    chromosome = chromosomes{j};
        
    % Iterate over all resolutions
    for l = 1:numel(resolution_values)
        resolution = resolution_values{l};
        label = resolution_labels{l};
        
        % Construct file paths    
        tensor_file_path = sprintf('%s/hic_methy_%s_all_cells_tensors/%s_all_cells_tensor.h5', output_directory, label, chromosome);
        output_file_U = sprintf('%s/tensor_%s_cell_type_single_modality_factors/%s_genomic_factors_modality0.h5', output_directory, label, chromosome);
        output_file_V = sprintf('%s/tensor_%s_cell_type_single_modality_factors/%s_sample_factors_modality1.h5', output_directory, label, chromosome);

        % Ensure output directories exist
        [output_dir_U, ~] = fileparts(output_file_U);
        [output_dir_V, ~] = fileparts(output_file_V);

        if ~exist(output_dir_U, 'dir')
            mkdir(output_dir_U);
        end

        if ~exist(output_dir_V, 'dir')
            mkdir(output_dir_V);
        end

        % Check if the tensor file exists before trying to load it
        if ~exist(tensor_file_path, 'file')
            fprintf('Tensor file %s does not exist. Skipping.\n', tensor_file_path);
            continue;
        end

        try
            % Load the tensor from the HDF5 file
            tensor = h5read(tensor_file_path, '/Combined_Tensor');
            tensor_size = size(tensor);
            fprintf('Processing tensor of size %s for %s.\n', mat2str(tensor_size), tensor_file_path);

            % Process each modality separately
            for modality = 1:2
                % Extract the modality-specific tensor and reshape
                modality_tensor = squeeze(tensor(:,modality,:,:)); % size n x n x 129
                modality_label = sprintf('%s_modality%d', label, modality);

                % Perform the tensor computation for the single modality
                model = struct;
                model.variables.u = randn(size(modality_tensor, 3), 4);
                model.variables.v = randn(size(modality_tensor, 1), 4);
                model.factors.U = {'u', @struct_nonneg};
                model.factors.V = {'v', @struct_nonneg};
                model.factorizations.myfac.data = modality_tensor;
                model.factorizations.myfac.cpd = {'V', 'U', 'U'};
                options.Display = 100;
                options.MaxIter = iterations;
                options.TolFun = 1e-15;
                options.TolX = 1e-15;
                sol = sdf_nls(model, options);

                % Extract the required fields from the nested structure
                U = sol.factors.U;
                V = sol.factors.V;

                % Create and write to the HDF5 files for this modality
                output_dataset_U = '/genomic_weights';
                output_dataset_V = '/sample_weights';

                output_file_U_mod = sprintf('%s/tensor_%s_cell_type_single_modality_factors/%s_genomic_factors_modality%d.h5', output_directory, label, chromosome, modality);
                output_file_V_mod = sprintf('%s/tensor_%s_cell_type_single_modality_factors/%s_sample_factors_modality%d.h5', output_directory, label, chromosome, modality);

                h5create(output_file_U_mod, output_dataset_U, size(U));
                h5write(output_file_U_mod, output_dataset_U, U);

                h5create(output_file_V_mod, output_dataset_V, size(V));
                h5write(output_file_V_mod, output_dataset_V, V);

                fprintf('Computed and saved results for modality %d to %s and %s.\n', modality, output_file_U_mod, output_file_V_mod);
            end
        catch ME
            fprintf('Error processing tensor for %s: %s\n', tensor_file_path, ME.message);
        end
    end
end

