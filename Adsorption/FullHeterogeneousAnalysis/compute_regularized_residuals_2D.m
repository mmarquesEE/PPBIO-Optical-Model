function residuals_aug = compute_regularized_residuals_2D(log_params, L, lambda, exp_settings, exp_data, model_config)
    % Calculates an an AUGMENTED residual vector that includes both the
    % data-fitting term and the regularization penalty term.
    
    p_linear = 10.^log_params;
    n_exp = length(exp_data);
    
    % Calculate the model-data residuals (this part is the same as before)
    model_signals_all_exp = cell(n_exp, 1);
    parfor exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [~, s_model] = run_single_experiment_2D_model(p_linear, setting, model_config);
        model_signals_all_exp{exp_idx} = s_model;
    end
    residuals_all_exp = cell(n_exp, 1);
    for exp_idx = 1:n_exp
        s_data = exp_data{exp_idx}.signals;
        s_model = model_signals_all_exp{exp_idx};
        residuals_all_exp{exp_idx} = s_model(:) - s_data(:);
    end
    model_residuals = cell2mat(residuals_all_exp);
    
    % --- NEW: Calculate the regularization penalty term ---
    % This penalizes solutions where adjacent parameters are very different.
    reg_penalty = lambda * L * p_linear;
    
    % --- NEW: Create the augmented residual vector ---
    residuals_aug = [model_residuals; reg_penalty];
end