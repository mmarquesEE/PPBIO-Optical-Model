function residuals = compute_residuals_2D_model(log_params, exp_settings, exp_data, model_config)
    % This function calculates the difference between the model output (using
    % the 2D parameterization) and the "experimental" data.
    
    p_linear = 10.^log_params;
    n_exp = length(exp_data);
    
    % Pre-allocate a cell array for model outputs
    model_signals_all_exp = cell(n_exp, 1);
    
    % Using parfor here can significantly speed up the optimization if you
    % have multiple experiments to fit simultaneously.
    parfor exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [~, s_model] = run_single_experiment_2D_model(p_linear, setting, model_config);
        model_signals_all_exp{exp_idx} = s_model;
    end
    
    % --- Combine residuals from all experiments ---
    residuals_all_exp = cell(n_exp, 1);
    for exp_idx = 1:n_exp
        s_data = exp_data{exp_idx}.signals;
        s_model = model_signals_all_exp{exp_idx};
        residuals_all_exp{exp_idx} = s_model(:) - s_data(:);
    end
    
    residuals = cell2mat(residuals_all_exp);
end