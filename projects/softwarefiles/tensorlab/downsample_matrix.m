% Define the downsample_matrix function
function downsampled_matrix = downsample_matrix(matrix, percent)
    [rows, cols] = size(matrix);
    num_entries = rows * cols;
    num_selected_entries = floor(num_entries * percent / 100);

    % Flatten the matrix and randomly choose the selected entries
    flattened_matrix = matrix(:);
    indices = randperm(num_entries, num_selected_entries);
    downsampled_matrix = zeros(num_entries, 1);
    downsampled_matrix(indices) = flattened_matrix(indices);

    % Reshape the downsampled matrix back to the original shape
    downsampled_matrix = reshape(downsampled_matrix, rows, cols);
end
