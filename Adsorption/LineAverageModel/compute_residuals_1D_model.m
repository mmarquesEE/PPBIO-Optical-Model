function all_residuals = compute_residuals_1D_model(log_params, exp_settings, exp_data, model_config)
% Computes residuals and includes robust error handling for simulation failures.
% This version is optimized to preallocate the results vector for performance.

    p_1D = 10.^log_params;

    % --- START: PREALLOCATION LOGIC ---

    % 1. First, calculate the total number of residual points across all experiments.
    total_num_residuals = 0;
    for exp_idx = 1:length(exp_settings)
        % numel() gets the total number of elements in the data matrix
        total_num_residuals = total_num_residuals + numel(exp_data{exp_idx}.signals);
    end
    
    % 2. Preallocate the full residuals vector with zeros.
    all_residuals = zeros(total_num_residuals, 1);
    
    % 3. Initialize an index to track our position in all_residuals.
    current_idx_start = 1;

    % --- END: PREALLOCATION LOGIC ---


    % --- Part 2: Model-Data Residuals with Error Handling ---
    for exp_idx = 1:length(exp_settings)
        setting = exp_settings(exp_idx);
        s_data_matrix = exp_data{exp_idx}.signals;
        
        try
            [~, s_sim_matrix] = run_single_experiment_1D_model(p_1D, setting, model_config);
            if ~isequal(size(s_sim_matrix), size(s_data_matrix))
                error('Simulation output size mismatch.');
            end
            res_matrix = s_sim_matrix - s_data_matrix;
        catch ME
            fprintf('Warning: Simulation failed. Penalizing parameter set. Error: %s\n', ME.message);
            
            penalty_value = 1e6 * (1 + norm(s_data_matrix(:)));
            res_matrix = penalty_value * ones(size(s_data_matrix));
        end
        
        % --- START: MODIFIED RESULT HANDLING ---
        
        % Calculate the number of elements for this specific experiment
        num_res_in_exp = numel(res_matrix);
        
        % Define the range in the preallocated vector to fill
        idx_range = current_idx_start : (current_idx_start + num_res_in_exp - 1);
        
        % Place the current residuals into the correct slice, ensuring it's a column
        all_residuals(idx_range) = res_matrix(:);
        
        % Update the starting index for the next iteration
        current_idx_start = current_idx_start + num_res_in_exp;
        
        % --- END: MODIFIED RESULT HANDLING ---
    end
end