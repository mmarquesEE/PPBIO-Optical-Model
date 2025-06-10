% Main script to run the Tikhonov regularization scenarios
clearvars; close all; clc;

fprintf('This script will demonstrate the necessity of Tikhonov regularization for parameter estimation.\n');
fprintf('It will run three scenarios and generate plots for each.\n\n');

% --- SCENARIO 1: Ideal case ---
% No noise, no regularization. We expect a perfect parameter recovery.
fprintf('--- RUNNING SCENARIO 3: Noise, with Regularization ---\n');
run_estimation_scenario('3: Noise, Regularization', 0.02, 5e-3);
fprintf('--- SCENARIO 3 COMPLETE ---\n\n');


function run_estimation_scenario(scenario_title, noise_level, lambda_tikhonov)
    % Shared parameters
    gridN_x = 10; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [5,5]; ads_y_range = [2,3];
    ads_nx = ads_x_range(2) - ads_x_range(1) + 1;
    ads_ny = ads_y_range(2) - ads_y_range(1) + 1;
    num_ads_cells = ads_nx * ads_ny;
    
    % Homogeneous parameters
    homog_params = [9.4e3, 0.0078, 1.0]; % kon, koff, smax_total
    smax_per_cell = homog_params(3) / num_ads_cells;
    
    % Fixed parameters
    D_coeff = 6e-3;
    ru_to_m = 1e-6;
    base_max_velocity = 8.3;
    
    % Example pulse parameters
    n_exp = 1;
    % Base parameters
    base_max_velocity = 8.3;
    T1 = 2000; T2 = 4000; T3 = 3*T1;
    t_total = 4*T1;
    c_diss = 0;
    c1 = 3.3e-6;
    c2 = 1e-6;

    % Generate 12 orthogonal experiments
    exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);

    % Display generated experiments
    for i = 1:length(exp_settings)
        fprintf('\nExperiment %d:\n', i);
        fprintf('  Concentrations: [%.2e, %.2e, %.2e]\n', exp_settings(i).pulse_concs);
        fprintf('  Velocity: %.2f\n', exp_settings(i).max_velocity);
    end
    
    % Create heterogeneous parameters
    [kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
        create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
    % Preallocate for identifiability analysis
    s_obs_final_base = zeros(n_exp,1); % Final s_obs values for baseline (each experiment)
    
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        fprintf('\nRunning Experiment %d\n', exp_idx);
        
        % Create homogeneous grids
        [kon_grid_homog, koff_grid_homog, smax_grid_homog] = ...
            create_homogeneous_grids(homog_params, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
        
        % Create velocity profile
        velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);
        
        % Time breaks and concentrations
        t_breaks = [0, setting.pulse_times, setting.t_total];
        concentrations = [setting.pulse_concs, setting.c_diss];
        
        % Initial condition for s (empty)
        s0_grid = zeros(gridN_x, gridN_y, gridN_z);
        
        % Homogeneous simulation
        [t_homog, ~, s_homog, K_homog, Q_homog] = simulate_3d_flow_model_with_pulses(...
            gridN_x, gridN_y, gridN_z, kon_grid_homog, koff_grid_homog, smax_grid_homog, ...
            velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid);
        
        % Heterogeneous simulation
        [t_heterog, ~, s_heterog, K_heterog, Q_heterog] = simulate_3d_flow_model_with_pulses(...
            gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
            velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid);
        
        % ============== KERNEL-BASED COMPOSITE BEHAVIOR ================
        % Extract adsorption region parameters
        kon_homog_ads = kon_grid_homog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        smax_homog_ads = smax_grid_homog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        kon_heterog_ads = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        smax_heterog_ads = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        
        % Extract adsorption region for state variables
        K_homog_ads = K_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        Q_homog_ads = Q_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        s_homog_ads = s_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        
        K_heterog_ads = K_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        Q_heterog_ads = Q_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        s_heterog_ads = s_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        
        % Get initial surface concentration (t=0)
        s0_homog_ads = squeeze(s_homog_ads(1, :, :, :));
        s0_heterog_ads = squeeze(s_heterog_ads(1, :, :, :));
        
        % Compute alpha (kon * smax) for each cell
        alpha_homog_ads = kon_homog_ads .* smax_homog_ads;
        alpha_heterog_ads = kon_heterog_ads .* smax_heterog_ads;
        
        % Reshape for broadcasting
        alpha_homog_ads = reshape(alpha_homog_ads, [1, size(alpha_homog_ads)]);
        alpha_heterog_ads = reshape(alpha_heterog_ads, [1, size(alpha_heterog_ads)]);
        s0_homog_ads = reshape(s0_homog_ads, [1, size(s0_homog_ads)]);
        s0_heterog_ads = reshape(s0_heterog_ads, [1, size(s0_heterog_ads)]);
        
        % Compute decay term: s0 * exp(-Q)
        decay_term_homog = s0_homog_ads .* exp(-Q_homog_ads);
        decay_term_heterog = s0_heterog_ads .* exp(-Q_heterog_ads);
        
        % Compute kernel term: alpha * K
        kernel_term_homog = alpha_homog_ads .* K_homog_ads;
        kernel_term_heterog = alpha_heterog_ads .* K_heterog_ads;
        
        % Sum over adsorption region
        decay_sum_homog = squeeze(sum(decay_term_homog, [2,3,4]));
        kernel_sum_homog = squeeze(sum(kernel_term_homog, [2,3,4]));
        s_obs_homog = decay_sum_homog + kernel_sum_homog;
        
        decay_sum_heterog = squeeze(sum(decay_term_heterog, [2,3,4]));
        kernel_sum_heterog = squeeze(sum(kernel_term_heterog, [2,3,4]));
        s_obs_heterog = decay_sum_heterog + kernel_sum_heterog;
        
        % Direct s_obs for verification
        s_obs_homog_direct = compute_s_obs(s_homog, ads_x_range, ads_y_range, ads_layer);
        s_obs_heterog_direct = compute_s_obs(s_heterog, ads_x_range, ads_y_range, ads_layer);
        
        % Store final s_obs value for identifiability analysis
        s_obs_final_base(exp_idx) = s_obs_heterog_direct(end);
        
        % Calculate discrepancy
        discrepancy = trapz(t_heterog, (s_obs_homog_direct - s_obs_heterog_direct).^2);
        fprintf('Composite behavior discrepancy: %.2e\n', discrepancy);
        
        % ====================== PLOTTING =========================
        % Plot composite behavior
        figure('Position', [100, 100, 1200, 800]);
        
        % Composite behavior comparison
        subplot(2,2,1);
        plot(t_homog, s_obs_homog_direct, 'b-', 'LineWidth', 2); hold on;
        plot(t_heterog, s_obs_heterog_direct, 'r--', 'LineWidth', 1.5);
        y_lims = [min([s_obs_homog_direct; s_obs_heterog_direct]), max([s_obs_homog_direct; s_obs_heterog_direct])];
        for i = 1:length(setting.pulse_times)
            line([setting.pulse_times(i), setting.pulse_times(i)], y_lims, ...
                'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        end
        xlabel('Time (s)'); ylabel('s_{obs}(t)');
        title(sprintf('Exp %d: Composite Behavior (v_{max}=%.1f)', exp_idx, setting.max_velocity));
        legend('Homogeneous', 'Heterogeneous (5% var)', 'Location', 'best');
        grid on;
        
        % Homogeneous decomposition
        subplot(2,2,3);
        plot(t_homog, s_obs_homog_direct, 'k-', 'LineWidth', 2); hold on;
        plot(t_homog, decay_sum_homog, 'b--', 'LineWidth', 1.5);
        plot(t_homog, kernel_sum_homog, 'r--', 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel('s_{obs}');
        legend('Total', 'Decay Term', 'Kernel Term');
        title('Homogeneous Case: Signal Decomposition');
        grid on;
        
        % Heterogeneous decomposition
        subplot(2,2,4);
        plot(t_heterog, s_obs_heterog_direct, 'k-', 'LineWidth', 2); hold on;
        plot(t_heterog, decay_sum_heterog, 'b--', 'LineWidth', 1.5);
        plot(t_heterog, kernel_sum_heterog, 'r--', 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel('s_{obs}');
        legend('Total', 'Decay Term', 'Kernel Term');
        title('Heterogeneous Case: Signal Decomposition');
        grid on;
        
        % Discrepancy plot
        subplot(2,2,2);
        plot(t_heterog, s_obs_homog_direct - s_obs_heterog_direct, 'm-', 'LineWidth', 1.5);
        xlabel('Time (s)'); ylabel('Difference');
        title(sprintf('Homog - Heterog\nDiscrepancy: %.2e', discrepancy));
        grid on;
        
        set(gcf, 'Name', sprintf('Experiment %d Results', exp_idx));
    end
    
    % ============== IDENTIFIABILITY ANALYSIS ================
    fprintf('\nStarting Identifiability Analysis...\n');
    
    % Extract adsorption region parameters
    kon_ads = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    koff_ads = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    smax_ads = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    % Flatten parameters into vector: [kon; koff; smax]
    p_true = [kon_ads(:); koff_ads(:); smax_ads(:)]; 
    N_params = length(p_true);
    
    % Preallocate combined Jacobian
    J_combined = [];
    
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        fprintf('Computing Jacobian for Experiment %d...\n', exp_idx);
        
        % Get time vector and composite signal for this experiment
        [~, s_obs] = run_single_experiment(...
            gridN_x, gridN_y, gridN_z, kon_ads, koff_ads, smax_ads,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        
        % Compute Jacobian via finite differences
        J_exp = zeros(length(s_obs), N_params);
        h = 1e-10;  % A good, robust choice for the relative step factor
        
        for k = 1:N_params
            p_pert = p_true;
            
            p_pert(k) = p_true(k) + h;
        
            % Split parameters... (rest of your code is the same)
            kon_pert = reshape(p_pert(1:numel(kon_ads)), size(kon_ads));
            koff_pert = reshape(p_pert(numel(kon_ads)+1:numel(kon_ads)+numel(koff_ads)), size(koff_ads));
            smax_pert = reshape(p_pert(numel(kon_ads)+numel(koff_ads)+1:end), size(smax_ads));
            
            % Run simulation...
            [~, s_pert] = run_single_experiment(...
                gridN_x, gridN_y, gridN_z, kon_pert, koff_pert, smax_pert,...
                ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
            
            % Compute derivative using the actual step_size used
            J_exp(:, k) = (s_pert - s_obs) / h;
        end
        
        % Stack Jacobians vertically
        J_combined = [J_combined; J_exp];
    end
    
    % ======== IDENTIFIABILITY CONDITION VERIFICATION ========
    % Compute SVD
    [U, S, V] = svd(J_combined, 'econ');
    svals = diag(S);
    
    % Calculate metrics
    rankJ = sum(svals > 1e-6 * max(svals));
    cond_number = cond(J_combined);
    identifiability_ratio = rankJ / N_params;
    
    fprintf('\nIdentifiability Analysis Results:\n');
    fprintf('Total Parameters: %d\n', N_params);
    fprintf('Rank of Combined Jacobian: %d\n', rankJ);
    fprintf('Condition Number: %.2e\n', cond_number);
    fprintf('Identifiability Ratio: %.2f\n', identifiability_ratio);
    
    if rankJ == N_params
        fprintf('--> SUFFICIENT RANK: Exact identifiability possible\n');
    else
        fprintf('--> INSUFFICIENT RANK: Identifiability not guaranteed\n');
    end
    
    % ====================== PLOTTING ========================
    % Singular value plot
    figure('Position', [100, 100, 1200, 500]);
    subplot(1,2,1);
    semilogy(svals, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    yline(1e-6 * max(svals), 'r--', 'Threshold', 'LineWidth', 1.5);
    xlabel('Singular Value Index');
    ylabel('Magnitude (log scale)');
    title('Singular Values of Combined Jacobian');
    grid on;
    legend('Singular Values', 'Rank Threshold');
    
    % Parameter sensitivity plot
    param_sensitivity = sqrt(sum(V.^2, 1));
    [~, param_order] = sort(param_sensitivity, 'descend');
    
    subplot(1,2,2);
    bar(param_sensitivity(param_order));
    xlabel('Parameter Index (sorted)');
    ylabel('Sensitivity Norm');
    title('Parameter Sensitivity Ranking');
    grid on;
    % --- SAVE FIGURE ---
    fig_filename = sprintf('Identifiability_Analysis_%s.png', strrep(strrep(scenario_title, ':', ''), ' ', '_'));
    saveas(gcf, fig_filename);
    fprintf('Saved figure to %s\n', fig_filename);
    % Add text annotations
    annotation('textbox', [0.15, 0.15, 0.3, 0.1], 'String', ...
        sprintf('Rank = %d/%d\nCond = %.1e', rankJ, N_params, cond_number), ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
    
    % ============== ENHANCED PARAMETER IDENTIFICATION ================    
    % Extract true heterogeneous parameters
    true_kon = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    true_koff = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    true_smax = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    % Generate "experimental" data from true parameters
    exp_data = cell(n_exp,1);
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [t, s_obs] = run_single_experiment(...
            gridN_x, gridN_y, gridN_z, true_kon, true_koff, true_smax,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        exp_data{exp_idx} = s_obs.* (1 + noise_level*randn(size(s_obs)));
    end
    
    % Prepare logarithmic parameterization
    num_ads_cells = numel(true_kon);
    homog_kon = homog_params(1);
    homog_koff = homog_params(2);
    homog_smax_cell = homog_params(3) / num_ads_cells;
    
    % Initial guess in log space
    rng(2); % For reproducible initialization
    p0_kon = log10(homog_kon) + 1*randn(num_ads_cells,1);
    p0_koff = log10(homog_koff) + 1*randn(num_ads_cells,1);
    p0_smax = log10(homog_smax_cell) + 1*randn(num_ads_cells,1);
    p0 = [p0_kon; p0_koff; p0_smax];
    
    % Bounds in log scale
    lb_kon = log10(1e2); ub_kon = log10(1e5); % kon between 100 and 100,000
    lb_koff = log10(1e-4); ub_koff = log10(1e0); % koff between 0.0001 and 1
    lb_smax = log10(1e-3); ub_smax = log10(1e1); % smax per cell between 0.001 and 10
    lb = [repmat(lb_kon, num_ads_cells, 1); ...
          repmat(lb_koff, num_ads_cells, 1); ...
          repmat(lb_smax, num_ads_cells, 1)];
    ub = [repmat(ub_kon, num_ads_cells, 1); ...
          repmat(ub_koff, num_ads_cells, 1); ...
          repmat(ub_smax, num_ads_cells, 1)];
    
   % --- Tikhonov Regularization Setup ---
    if lambda_tikhonov > 0
        fprintf('Using SPATIAL regularization (lambda = %.2e)\n', lambda_tikhonov);
        % ** NEW: Create the spatial regularization operator L **
        L_matrix_reg = create_spatial_regularization_operator(ads_nx, ads_ny);
        % For pure smoothness, the reference vector is zero.
        log_params_ref_for_reg = zeros(N_params, 1);
    else
        fprintf('No regularization being used (lambda = 0)\n');
        % L matrix is not needed, but must be passed to the function
        L_matrix_reg = []; 
        log_params_ref_for_reg = [];
    end

    % --- Create a handle to the plotter function with all necessary data ---
    model_config.gridN_x = gridN_x; model_config.gridN_y = gridN_y; model_config.gridN_z = gridN_z;
    model_config.ads_x_range = ads_x_range; model_config.ads_y_range = ads_y_range; model_config.ads_layer = ads_layer;
    model_config.D_coeff = D_coeff; model_config.ru_to_m = ru_to_m;
    % Ensure this variable is available for error calculation
    true_kon_ads_region_for_error = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    true_koff_ads_region_for_error = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    true_smax_ads_region_for_error = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);

    true_params_vec = [true_kon_ads_region_for_error(:); true_koff_ads_region_for_error(:); true_smax_ads_region_for_error(:)]; 
    % --- Setup for Video Recording ---
    video_filename = sprintf('optimization_animation+%s.mp4', scenario_title);
    try
        video_obj = VideoWriter(video_filename, 'MPEG-4');
        video_obj.FrameRate = 5;  % Adjust FrameRate for desired speed (e.g., 5-10)
        video_obj.Quality = 95; % Adjust Quality (0-100, higher is better)
    catch ME
        warning('VideoWriter could not be created. Video will not be saved. Error: %s', ME.message);
        video_obj = []; % Set to empty if creation fails
    end
    optim_plot_fun = @(log_params, optimValues, state) optimPlotter(log_params, optimValues, state, ...
                                                              true_params_vec, ...
                                                              exp_data, ...
                                                              exp_settings, ...
                                                              model_config, ...
                                                              size(true_kon), size(true_koff), size(true_smax), ...
                                                              video_obj); % Pass the video object as the last argument

    % --- Add the 'OutputFcn' to your optimization options ---
    optim_opts = optimoptions('lsqnonlin', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', ...
        'MaxIterations', 150, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'MaxFunctionEvaluations',3000, ...
        'UseParallel',true,...
        'FiniteDifferenceType', 'central', ...
        'OutputFcn', optim_plot_fun); % This now points to our video-enabled plotter
    
    % Residual function with Tikhonov regularization
    % Ensure n_exp is correctly passed (it's the number of experiments used for fitting)
    num_fitting_experiments = n_exp; % Assuming you use all 'n_exp' for fitting

    residual_fun_tikhonov = @(log_params) compute_residuals_log_tikhonov(...
        log_params, exp_settings(1:num_fitting_experiments), exp_data(1:num_fitting_experiments), ... % Use only the selected/fitted experiments
        gridN_x, gridN_y, gridN_z, ...
        ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m, ...
        num_fitting_experiments, ... % Pass n_exp correctly
        lambda_tikhonov, log_params_ref_for_reg, L_matrix_reg); % Pass regularization params
    
    fprintf('\nStarting Enhanced Parameter Identification with Tikhonov Regularization (lambda=%.1e)...\n', lambda_tikhonov);
    
    % Run optimization in log space
    [opt_log_params_reg, resnorm_reg, residual_reg, exitflag_reg] = lsqnonlin(...
    residual_fun_tikhonov, p0, lb, ub, optim_opts);
    
    % Convert optimized parameters back to linear space
    opt_params_reg = 10.^opt_log_params_reg;
    
    % Split into parameter groups
    opt_kon_reg = opt_params_reg(1:num_ads_cells);
    opt_koff_reg = opt_params_reg(num_ads_cells+1:2*num_ads_cells);
    opt_smax_reg = opt_params_reg(2*num_ads_cells+1:end);
    
    % Reshape to match adsorption region
    opt_kon_reg = reshape(opt_kon_reg, size(true_kon));
    opt_koff_reg = reshape(opt_koff_reg, size(true_koff));
    opt_smax_reg = reshape(opt_smax_reg, size(true_smax));
    
    % True parameters as vector for comparison (already have true_params_vec)
    
    % Analyze regularized results
    fprintf('\nOptimization Results (Tikhonov Regularized):\n');
    fprintf('Final Residual Norm (augmented): %.4e\n', norm(residual_reg)); % This includes reg term
    % To get model-only residual norm:
    num_model_residuals = length(residual_reg) - N_params; % If L=I
    if lambda_tikhonov > 0 && size(L_matrix_reg,1) == N_params % Assuming L_matrix_reg results in N_params rows for reg term
        model_only_resnorm = norm(residual_reg(1:num_model_residuals));
        fprintf('Model-Data Residual Norm (unaugmented part): %.4e\n', model_only_resnorm);
    end
    fprintf('Exit Flag: %d\n', exitflag_reg);
        
    
    % Calculate parameter errors for regularized solution
    param_errors_reg = abs(opt_params_reg - true_params_vec) ./ true_params_vec;
    fprintf('\nParameter Recovery Accuracy (Tikhonov Regularized):\n');
    fprintf('Mean Relative Error: %.2f%%\n', 100*mean(param_errors_reg(~isinf(param_errors_reg) & ~isnan(param_errors_reg))));
    fprintf('Max Relative Error: %.2f%%\n', 100*max(param_errors_reg(~isinf(param_errors_reg) & ~isnan(param_errors_reg))));
    init_params = 10.^p0; % Initial guess in linear space
    % Plot parameter recovery for regularized solution
    plot_parameter_recovery(true_params_vec, opt_params_reg, init_params, ... % init_params is 10.^p0
        size(true_kon), size(true_koff), size(true_smax),scenario_title);
    sgtitle('Parameter Recovery Results (Tikhonov Regularized)'); % Add to distinguish plot

    % Plot predicted vs "experimental" signals using regularized parameters
    plot_signal_predictions([opt_kon_reg(:); opt_koff_reg(:); opt_smax_reg(:)], ...
        exp_settings(1:num_fitting_experiments), exp_data(1:num_fitting_experiments), ...
        gridN_x, gridN_y, gridN_z, ...
        ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m, num_fitting_experiments,scenario_title);
    sgtitle('Signal Prediction vs Experimental Data (Tikhonov Regularized Parameters)');
end

function exp_settings = generate_experiments(M, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2)
    % Generates M experiments with orthogonal concentration profiles and flow velocities
    %
    % Inputs:
    %   M - Number of experiments
    %   base_max_velocity - Reference flow velocity (e.g., 8.3)
    %   T1, T2, T3 - Fixed pulse times
    %   t_total - Total experiment duration
    %   c_diss - Dissociation concentration
    %   c1, c2 - Base concentrations
    
    exp_settings = struct(...
        'pulse_times', {}, ...
        'pulse_concs', {}, ...
        't_total', {}, ...
        'max_velocity', {}, ...
        'c_diss', {} ...
    );
    
    % Generate orthogonal concentration pairs using polar coordinates
    angles = linspace(0, pi/2, M);  % Cover quadrant for positive concentrations
    factors = linspace(0.3, 3, M);  % Concentration scaling factors
    
    for i = 1:M
        if M ==1
            M=2;
        end
        % Create orthogonal concentration profiles
        conc_factor1 = factors(ceil(i/2)) * cos(angles(i));
        conc_factor2 = factors(ceil(i/2)) * sin(angles(i));
        
        % Ensure minimum concentration variation
        min_conc = 0.1 * min(c1, c2);
        conc1 = max(c1 * (0.5 + conc_factor1), min_conc);
        conc3 = max(c2 * (0.5 + conc_factor2), min_conc);
        
        % Create velocity profile (logarithmic spacing)
        vel_min = 0.1 * base_max_velocity;
        vel_max = 5.0 * base_max_velocity;
        velocity = exp(log(vel_min) + (i-1)/(M-1) * (log(vel_max) - log(vel_min)));
        
        % Special patterns for every 3rd experiment
        if mod(i,3) == 0
            exp_settings(i).pulse_concs = [conc1, c2, conc3];  % Middle pulse active
        elseif mod(i,4) == 0
            exp_settings(i).pulse_concs = [c1, 0, conc3];      % First pulse fixed
        else
            exp_settings(i).pulse_concs = [conc1, 0, conc3];   % Standard pattern
        end
        
        % Assign common parameters
        exp_settings(i).pulse_times = [T1, T2, T3];
        exp_settings(i).t_total = t_total;
        exp_settings(i).max_velocity = velocity;
        exp_settings(i).c_diss = c_diss;
    end
end
% ================== NEW HELPER FUNCTIONS ==================
% --- Modified compute_residuals_log function ---
function L_full = create_spatial_regularization_operator(nx, ny)
    % Creates a sparse finite difference operator L = [Dx; Dy] for a single
    % parameter field, then combines them for all 3 parameter types.
    % Assumes column-major ordering of parameters: p = [p11, p21, ..., pNx1, p12, ...]
    
    N = nx * ny;
    
    % X-derivative (differences between i+1 and i)
    e_x = ones(N, 1);
    Dx = spdiags([-e_x, e_x], [0, 1], N, N);
    % Remove wrap-around connections for each column
    for i = 1:(ny-1)
        Dx(i*nx, i*nx + 1) = 0;
    end
    Dx(ny*nx, :) = []; % Remove last row which is all zeros

    % Y-derivative (differences between j+1 and j)
    e_y = ones(N, 1);
    Dy = spdiags([-e_y, e_y], [0, nx], N, N);
    % Remove connections for the last ny-nx rows
    Dy(N-nx+1:end, :) = [];
    
    L_space = [Dx; Dy];

    % Combine for all 3 parameter types (kon, koff, smax)
    L_full = blkdiag(L_space, L_space, L_space);
    
    fprintf('Created spatial regularization operator L of size %d x %d\n', ...
        size(L_full, 1), size(L_full, 2));
end


function residuals_aug = compute_residuals_log_tikhonov(...
    log_params, exp_settings, exp_data, ...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m, n_exp, ...
    lambda_reg, log_params_ref, L_matrix)
    
    params = 10.^log_params;
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    ads_param_shape = [ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1];
    
    kon_ads = reshape(params(1:num_ads_cells), ads_param_shape);
    koff_ads = reshape(params(num_ads_cells+1:2*num_ads_cells), ads_param_shape);
    smax_ads = reshape(params(2*num_ads_cells+1:end), ads_param_shape);
    
    % --- Model-Data Mismatch (Fidelity Term) ---
    model_residuals = [];
    sim_data = cell(n_exp, 1);

    parfor i = 1:n_exp
        setting = exp_settings(i);
        [~, s_sim] = run_single_experiment(...
            nx, ny, nz, kon_ads, koff_ads, smax_ads,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        sim_data{i} = s_sim;
    end
    
    for i = 1:n_exp
        model_residuals = [model_residuals; (sim_data{i} - exp_data{i})];
    end

    % --- Tikhonov Regularization Term ---
    if lambda_reg > 0 && ~isempty(L_matrix)
        % ** MODIFIED: Penalize spatial differences using the L matrix **
        % The penalty is applied to the log-parameters for better scaling.
        % The reference vector log_params_ref should be zeros for pure smoothness.
        reg_term_vector = sqrt(lambda_reg) * (L_matrix * (log_params(:) - log_params_ref(:)));
        
        % Augment the residual vector
        residuals_aug = [model_residuals; reg_term_vector];
    else
        % No regularization
        residuals_aug = model_residuals;
    end
end


% ================== EXISTING HELPER FUNCTIONS ==================
function [t, s_obs] = run_single_experiment(...
    nx, ny, nz, kon_ads, koff_ads, smax_ads,...
    ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m)
    
    % Create full parameter grids
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_ads;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_ads;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_ads;
    
    % Create velocity profile
    velocity_profile = create_velocity_profile(nz, setting.max_velocity);
    
    % Time breaks and concentrations
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    
    % Initial condition
    s0_grid = zeros(nx, ny, nz);
    
    % Run simulation
    [t, ~, s] = simulate_3d_flow_model_with_pulses(...
        nx, ny, nz, kon_grid, koff_grid, smax_grid, ...
        velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid);
    
    % Compute observed signal
    s_obs = compute_s_obs(s, ads_x_range, ads_y_range, ads_layer);
end

function plot_parameter_recovery(true_params, opt_params, initial_guess, sz_kon, sz_koff, sz_smax,scenario_title)
    num_kon = prod(sz_kon);
    num_koff = prod(sz_koff);
    
    % Extract parameter groups
    true_kon = true_params(1:num_kon);
    true_koff = true_params(num_kon+1:num_kon+num_koff);
    true_smax = true_params(num_kon+num_koff+1:end);
    
    opt_kon = opt_params(1:num_kon);
    opt_koff = opt_params(num_kon+1:num_kon+num_koff);
    opt_smax = opt_params(num_kon+num_koff+1:end);
    
    init_kon = initial_guess(1:num_kon);
    init_koff = initial_guess(num_kon+1:num_kon+num_koff);
    init_smax = initial_guess(num_kon+num_koff+1:end);
    
    % Create figure
    figure('Position', [100, 100, 1200, 900]);
    
    % Plot kon recovery
    subplot(3,2,1);
    plot(true_kon, 'ro', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_kon, 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_kon, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Binding Rate (k_{on}) Recovery');
    legend('True', 'Initial Guess', 'Recovered');
    ylabel('Value');
    grid on;
    
    subplot(3,2,2);
    error_kon = abs(opt_kon - true_kon) ./ true_kon;
    bar(100*error_kon);
    title('Relative Error in k_{on}');
    ylabel('Error (%)');
    ylim([0, 50]);
    grid on;
    
    % Plot koff recovery
    subplot(3,2,3);
    plot(true_koff, 'ro', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_koff, 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_koff, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Unbinding Rate (k_{off}) Recovery');
    ylabel('Value');
    grid on;
    
    subplot(3,2,4);
    error_koff = abs(opt_koff - true_koff) ./ true_koff;
    bar(100*error_koff);
    title('Relative Error in k_{off}');
    ylabel('Error (%)');
    ylim([0, 50]);
    grid on;
    
    % Plot smax recovery
    subplot(3,2,5);
    plot(true_smax, 'ro', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_smax, 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_smax, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Maximum Binding (s_{max}) Recovery');
    xlabel('Parameter Index');
    ylabel('Value');
    grid on;
    
    subplot(3,2,6);
    error_smax = abs(opt_smax - true_smax) ./ true_smax;
    bar(100*error_smax);
    title('Relative Error in s_{max}');
    xlabel('Parameter Index');
    ylabel('Error (%)');
    ylim([0, 50]);
    grid on;
    
    sgtitle(['Parameter Recovery Results for Scenario: ' scenario_title], 'FontSize', 14, 'FontWeight', 'bold');
    % --- SAVE FIGURE ---
    fig_filename = sprintf('Parameter_Recovery_%s.png', strrep(strrep(scenario_title, ':', ''), ' ', '_'));
    saveas(gcf, fig_filename);
    fprintf('Saved figure to %s\n', fig_filename);
end

function plot_signal_predictions(opt_params, exp_settings, exp_data, ...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,n_exp,scenario_title)
    
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    kon_ads = reshape(opt_params(1:num_ads_cells), [ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1]);
    koff_ads = reshape(opt_params(num_ads_cells+1:2*num_ads_cells), size(kon_ads));
    smax_ads = reshape(opt_params(2*num_ads_cells+1:end), size(kon_ads));
    
    figure('Position', [100, 100, 1200, 400 * n_exp]);
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [t, s_sim] = run_single_experiment(...
            nx, ny, nz, kon_ads, koff_ads, smax_ads,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        
        subplot(3,1,exp_idx);
        plot(t, exp_data{exp_idx}, 'b-', 'LineWidth', 2); hold on;
        plot(t, s_sim, 'r--', 'LineWidth', 1.5);
        title(sprintf('Experiment %d: Signal Comparison (v_{max}=%.1f)', exp_idx, setting.max_velocity));
        xlabel('Time (s)');
        ylabel('s_{obs}(t)');
        legend('"Experimental"', 'Recovered Parameters', 'Location', 'best');
        grid on;
        
        % Add pulse indicators
        y_lims = ylim;
        for i = 1:length(setting.pulse_times)
            line([setting.pulse_times(i), setting.pulse_times(i)], y_lims, ...
                'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        end
    end
    sgtitle(['Signal Fit vs. Data for Scenario: ' scenario_title], 'FontSize', 14, 'FontWeight', 'bold');
    % --- SAVE FIGURE ---
    fig_filename = sprintf('Signal_Fit_%s.png', strrep(strrep(scenario_title, ':', ''), ' ', '_'));
    saveas(gcf, fig_filename);
    fprintf('Saved figure to %s\n', fig_filename);
end

function [kon_grid, koff_grid, smax_grid] = create_heterogeneous_grids_from_ads(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_ads, koff_ads, smax_ads)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_ads;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_ads;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_ads;
end

function s_obs = compute_s_obs(s_grid, ads_x_range, ads_y_range, ads_layer)
    ads_cells = s_grid(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    s_obs = squeeze(sum(ads_cells, [2,3,4]));
end

function [kon_grid, koff_grid, smax_grid] = create_homogeneous_grids(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    kon = params(1);
    koff = params(2);
    smax_total = params(3);
    smax_per_cell = smax_total / ((ads_x_range(2)-ads_x_range(1)+1)*(ads_y_range(2)-ads_y_range(1)+1));
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_per_cell;
end

function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz-1);
    h = nz-1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1,1,nz]);
end

function [kon_grid, koff_grid, smax_grid] = create_ground_truth_heterogeneity(nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    kon_base = 9.4e3;
    koff_base = 0.0078;
    smax_total = 1.0;
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    rng(42);
    kon_vals = kon_base * (1 + 0.05*randn(ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1));
    koff_vals = koff_base * (1 + 0.05*randn(size(kon_vals)));
    smax_vals = (smax_total/num_ads_cells) * (1 + 0.05*randn(size(kon_vals)));
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_vals;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_vals;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_vals;
end

function [t, c_s, s, K, Q] = simulate_3d_flow_model_with_pulses(...
    nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid)
    
    % Initialize state variables
    num_cells = nx * ny * nz;
    c_s0 = zeros(nx, ny, nz);
    c_s0(1, :, :) = concentrations(1); % Initial concentration
    s0 = s0_grid;
    Q0 = zeros(nx, ny, nz);
    R0 = zeros(nx, ny, nz);
    y0 = [c_s0(:); s0(:); Q0(:); R0(:)];
    
    % Setup ODE options
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
    
    % Preallocate results
    t_all = [];
    y_all = [];
    
    % Process each time segment
    num_segments = length(t_breaks) - 1;
    for seg = 1:num_segments
        t_start = t_breaks(seg);
        t_end = t_breaks(seg+1);
        c0_seg = concentrations(seg);
        
        % Determine time points for segment
        num_points = 1000;
        tspan = linspace(t_start, t_end, num_points);
        % Run simulation for segment
        [t_seg, y_seg] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, ...
            kon_grid, koff_grid, smax_grid, c0_seg, D_coeff, ru_to_m), tspan, y0, options);
        
        % Handle first segment specially
        if seg == 1
            t_all = t_seg;
            y_all = y_seg;
        else
            % Append results (skip first point to avoid duplicate)
            t_all = [t_all; t_seg(2:end)];
            y_all = [y_all; y_seg(2:end, :)];
        end
        
        % Update initial condition for next segment
        if seg < num_segments
            y0 = y_seg(end, :)';
            c_s_end = reshape(y0(1:num_cells), [nx, ny, nz]);
            s_end = reshape(y0(num_cells+1:2*num_cells), [nx, ny, nz]);
            c_s_end(1, :, :) = concentrations(seg+1);
            y0 = [c_s_end(:); s_end(:); y0(2*num_cells+1:end)]; % Preserve Q/R
        end
    end
    
    % Extract variables
    c_s = reshape(y_all(:, 1:num_cells), [length(t_all), nx, ny, nz]);
    s = reshape(y_all(:, num_cells+1:2*num_cells), [length(t_all), nx, ny, nz]);
    Q = reshape(y_all(:, 2*num_cells+1:3*num_cells), [length(t_all), nx, ny, nz]);
    R = reshape(y_all(:, 3*num_cells+1:end), [length(t_all), nx, ny, nz]);
    
    % Compute kernel
    K = exp(-Q) .* R;
    t = t_all;
end

function dydt = ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0, D_coeff, ru_to_m)
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);
    Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);
    R = reshape(y(3*num_cells + 1:4*num_cells), [nx, ny, nz]);
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);

    % Diffusion terms
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = (c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:));
    
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2);
    d2c_dz2(:,:,1) = c_s(:,:,2) - 2*c_s(:,:,1) + c_s(:,:,1);
    d2c_dz2(:,:,end) = c_s(:,:,end-1) - 2*c_s(:,:,end) + c_s(:,:,end-1);
    
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2);
    
    % Advection
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)));
    
    % Adsorption kinetics
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c0;
    dcsdt(1,:,:) = 0;
    
    % Compute dQ/dt and dR/dt
    dQdt = kon_grid .* c_s + koff_grid;
    dRdt = c_s .* exp(Q);

    % Combine all derivatives
    dydt = [dcsdt(:); dsdt(:); dQdt(:); dRdt(:)];
end
function stop = optimPlotter(log_params, optimValues, state, ...
                            true_params_vec, exp_data, exp_settings, model_config, ...
                            sz_kon, sz_koff, sz_smax, video_obj)
% optimPlotter is an output function for lsqnonlin to visualize optimization progress.
% It plots the evolution of parameters and sensorgram fit, and saves the animation to a video file.

    stop = false; % This function does not stop the optimization
    
    % Use persistent variables to store plot handles and static data
    persistent handles; 

    switch state
        case 'init'
            % On the first call, create the figure, axes, and open the video file
            fig = figure('Name', 'Optimization Progress', 'Position', [150, 150, 1600, 700]);
            
            % --- Setup Axes ---
            ax_kon = subplot(2, 3, 1);
            ax_koff = subplot(2, 3, 2);
            ax_smax = subplot(2, 3, 3);
            num_exp_to_plot = min(length(exp_settings), 3);
            ax_sigs = gobjects(1, num_exp_to_plot);
            for i = 1:num_exp_to_plot
                ax_sigs(i) = subplot(2, num_exp_to_plot, num_exp_to_plot + i);
            end
            
            % --- Store handles and static data ---
            handles.fig = fig;
            handles.ax_kon = ax_kon; handles.ax_koff = ax_koff; handles.ax_smax = ax_smax;
            handles.ax_sigs = ax_sigs;
            
            handles.true_params_vec = true_params_vec;
            handles.exp_data = exp_data;
            handles.exp_settings = exp_settings;
            handles.model_config = model_config;
            handles.sz_kon = sz_kon; handles.sz_koff = sz_koff; handles.sz_smax = sz_smax;
            handles.n_exp_to_plot = num_exp_to_plot;
            handles.video_obj = video_obj; % Store video object handle
            
            % --- Open video file for writing ---
            if ~isempty(handles.video_obj)
                try
                    open(handles.video_obj);
                catch ME
                    warning('Could not open video file for writing: %s', ME.message);
                    handles.video_obj = []; % Invalidate object if it fails
                end
            end
            
            % Call update function to draw and capture the initial frame
            updatePlotsAndVideo(log_params, optimValues, handles);

        case 'iter'
            % On each subsequent iteration, update the plots and save the frame
            if ishandle(handles.fig) % Check if the user has closed the figure
                updatePlotsAndVideo(log_params, optimValues, handles);
            else
                stop = true; % Stop optimization if figure is closed
                fprintf('Animation figure closed. Stopping optimization.\n');
            end
            
        case 'done'
            % Finalize and close the video file
            if isfield(handles, 'video_obj') && ~isempty(handles.video_obj)
                try
                    close(handles.video_obj);
                    fprintf('Video successfully saved to: %s\n', handles.video_obj.Filename);
                catch ME
                    warning('Could not finalize video file: %s', ME.message);
                end
            end
            if ishandle(handles.fig)
                sgtitle(handles.fig, 'Optimization Finished!', 'FontSize', 14, 'FontWeight', 'bold');
            end
    end

    % --- Nested function to handle plotting AND video frame writing ---
    function updatePlotsAndVideo(current_log_params, optimVals, h)
        
        current_params_linear = 10.^current_log_params;

        % --- Split parameters for plotting ---
        num_kon = prod(h.sz_kon);
        num_koff = prod(h.sz_koff);
        true_kon = h.true_params_vec(1:num_kon);
        true_koff = h.true_params_vec(num_kon+1:num_kon+num_koff);
        true_smax = h.true_params_vec(num_kon+num_koff+1:end);
        opt_kon = current_params_linear(1:num_kon);
        opt_koff = current_params_linear(num_kon+1:num_kon+num_koff);
        opt_smax = current_params_linear(num_kon+num_koff+1:end);

        % --- Update Parameter Plots ---
        cla(h.ax_kon); hold(h.ax_kon, 'on');
        plot(h.ax_kon, true_kon, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'True');
        plot(h.ax_kon, opt_kon, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Current');
        title(h.ax_kon, sprintf('k_{on} (Iter: %d)', optimVals.iteration));
        legend(h.ax_kon, 'Location', 'best'); grid(h.ax_kon, 'on'); hold(h.ax_kon, 'off');
        
        cla(h.ax_koff); hold(h.ax_koff, 'on');
        plot(h.ax_koff, true_koff, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(h.ax_koff, opt_koff, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
        title(h.ax_koff, sprintf('k_{off} (F-count: %d)', optimVals.funccount));
        grid(h.ax_koff, 'on'); hold(h.ax_koff, 'off');

        cla(h.ax_smax); hold(h.ax_smax, 'on');
        plot(h.ax_smax, true_smax, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        plot(h.ax_smax, opt_smax, 'g*', 'MarkerSize', 8, 'LineWidth', 1.5);
        title(h.ax_smax, sprintf('s_{max} (Residual: %.2e)', optimVals.resnorm));
        grid(h.ax_smax, 'on'); hold(h.ax_smax, 'off');

        % --- Update Sensorgram Plots ---
        kon_ads_current = reshape(opt_kon, h.sz_kon);
        koff_ads_current = reshape(opt_koff, h.sz_koff);
        smax_ads_current = reshape(opt_smax, h.sz_smax);
        num_model_residuals = numel([h.exp_data{:}]);
        model_residuals = optimVals.residual(1:num_model_residuals);
        
        % Reconstruct the full simulated signal for plotting
        s_sim_all = model_residuals + [h.exp_data{:}];
        start_idx = 1;
        for i_plot = 1:h.n_exp_to_plot
            num_pts = numel(h.exp_data{i_plot});
            s_sim_exp = s_sim_all(start_idx : start_idx + num_pts - 1);
            t = linspace(0, h.exp_settings(i_plot).t_total, num_pts); % Recreate time vector
            start_idx = start_idx + num_pts;
        
            ax = h.ax_sigs(i_plot);
            cla(ax); hold(ax, 'on');
            plot(ax, t, h.exp_data{i_plot}, 'b-', 'LineWidth', 2, 'DisplayName', 'Data');
            plot(ax, t, s_sim_exp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Current Fit');
            title(ax, sprintf('Sensorgram - Exp %d', i_plot));
            legend(ax, 'Location', 'best'); xlabel(ax, 'Time (s)'); ylabel(ax, 's_{obs}(t)');
            grid(ax, 'on'); hold(ax, 'off');
        end

        drawnow; % Force the figure window to update

        % --- Capture Frame for Video ---
        if ~isempty(h.video_obj)
            try
                frame = getframe(h.fig); % Capture the entire figure window
                writeVideo(h.video_obj, frame);
            catch ME
                warning('Could not write frame to video: %s', ME.message);
                % Invalidate object to stop trying on subsequent iterations
                h.video_obj = []; 
            end
        end
    end
end





% --- RUNNING SCENARIO 1: No Noise, No Regularization ---
% 
% Experiment 1:
%   Concentrations: [1.65e-06, 0.00e+00, 3.50e-06]
%   Velocity: 0.83
% 
% Running Experiment 1
% Composite behavior discrepancy: 4.58e+00
% 
% Starting Identifiability Analysis...
% Computing Jacobian for Experiment 1...
% 
% Identifiability Analysis Results:
% Total Parameters: 6
% Rank of Combined Jacobian: 6
% Condition Number: 9.39e+02
% Identifiability Ratio: 1.00
% --> SUFFICIENT RANK: Exact identifiability possible
% Saved figure to Identifiability_Analysis_1_No_Noise,_No_Regularization.png
% Warning: VideoWriter could not be created. Video will not be saved. Error: Unable to find file. Ensure file exists
% and path is valid. 
% > In run_identifiability_sweep>run_estimation_scenario (line 388)
% In run_identifiability_sweep (line 10) 
% 
% Starting Enhanced Parameter Identification with Tikhonov Regularization (lambda=0.0e+00)...
% 
%                                          Norm of      First-order 
%  Iteration  Func-count     f(x)          step          optimality
%      0         13         322.392                           466
%      1         26         23.1513         0.4787            489      
%      2         39        0.905494       0.130412           10.8      
%      3         52       0.0525238       0.132854           21.5      
%      4         65      0.00443217      0.0352917          0.135      
%      5         78      0.00134861       0.102728          0.245      
%      6         91     0.000372439      0.0522019           1.79      
%      7        104     0.000136328      0.0112515          0.614      
%      8        117     0.000136328      0.0165591          0.614      
%      9        130     9.33259e-05     0.00413978         0.0715      
%     10        143     7.76282e-05     0.00827956          0.507      
%     11        156     4.62979e-05     0.00827956          0.333      
%     12        169     3.33435e-05     0.00827956          0.379      
%     13        182      2.5121e-05     0.00835601          0.416      
%     14        195     1.05411e-05     0.00206989         0.0291      
%     15        208     8.06264e-06     0.00413978          0.103      
%     16        221     5.88096e-06     0.00413978          0.123      
%     17        234     4.26886e-06     0.00413978          0.134      
%     18        247      3.0857e-06     0.00389847          0.133      
%     19        260     1.49125e-06     0.00243828         0.0609      
%     20        273     1.49125e-06     0.00272703         0.0609      
%     21        286     1.08955e-06    0.000681758        0.00583      
%     22        299     9.63353e-07     0.00136352         0.0253      
%     23        312     8.29289e-07     0.00136352         0.0248      
%     24        325     7.38459e-07     0.00136352         0.0252      
%     25        338     6.72457e-07     0.00136352         0.0247      
%     26        351     6.20255e-07     0.00136352         0.0233      
%     27        364     6.16971e-07     0.00163291         0.0323      
%     28        377      5.1908e-07    0.000340879        0.00133      
%     29        390     5.04583e-07    0.000681758        0.00541      
%     30        403     5.04583e-07     0.00136352        0.00541      
%     31        416     4.93753e-07    0.000340879        0.00132      
%     32        429     4.80291e-07    0.000681758        0.00561      
%     33        442     4.80291e-07     0.00136352        0.00561      
%     34        455     4.69627e-07    0.000340879        0.00136      
%     35        468     4.56943e-07    0.000681758        0.00577      
%     36        481     4.56943e-07     0.00136352        0.00577      
%     37        494     4.46399e-07    0.000340879         0.0014      
%     38        507     4.34339e-07    0.000681758        0.00589      
%     39        520     4.34339e-07     0.00136352        0.00589      
%     40        533     4.23906e-07    0.000340879        0.00142      
%     41        546     4.12374e-07    0.000681758        0.00597      
%     42        559     4.12374e-07     0.00136352        0.00597      
%     43        572     4.02066e-07    0.000340879        0.00144      
%     44        585     3.91003e-07    0.000681758        0.00602      
%     45        598     3.91003e-07     0.00136352        0.00602      
%     46        611     3.80844e-07    0.000340879        0.00145      
%     47        624     3.70217e-07    0.000681758        0.00604      
%     48        637     3.70217e-07     0.00136352        0.00604      
%     49        650     3.60233e-07    0.000340879        0.00146      
%     50        663     3.50018e-07    0.000681758        0.00604      
%     51        676     3.50018e-07     0.00136352        0.00604      
%     52        689     3.40239e-07    0.000340879        0.00147      
%     53        702      3.3043e-07    0.000681758        0.00602      
%     54        715      3.3043e-07     0.00136352        0.00602      
%     55        728     3.20878e-07    0.000340879        0.00147      
%     56        741     3.11469e-07    0.000681758        0.00599      
%     57        754     3.11469e-07     0.00136352        0.00599      
%     58        767     3.02167e-07    0.000340879        0.00146      
%     59        780     2.93154e-07    0.000681758        0.00594      
%     60        793     2.81239e-07    0.000681758        0.00596      
%     61        806     2.81239e-07     0.00136352        0.00596      
%     62        819      2.7235e-07    0.000340879        0.00145      
%     63        832     2.63998e-07    0.000681758        0.00584      
%     64        845      2.5288e-07    0.000681758        0.00585      
%     65        858      2.5288e-07     0.00136352        0.00585      
%     66        871     2.44483e-07    0.000340879        0.00143      
%     67        884     2.36768e-07    0.000681758        0.00571      
%     68        897     2.26423e-07    0.000681758        0.00572      
%     69        910     2.26423e-07     0.00136352        0.00572      
%     70        923     2.18522e-07    0.000340879        0.00141      
%     71        936     2.11407e-07    0.000681758        0.00557      
%     72        949     2.01796e-07    0.000681758        0.00558      
%     73        962     2.01796e-07     0.00136352        0.00558      
%     74        975     1.94384e-07    0.000340879        0.00138      
%     75        988     1.87825e-07    0.000681758        0.00542      
%     76       1001       1.789e-07    0.000681758        0.00542      
%     77       1014       1.789e-07     0.00136352        0.00542      
%     78       1027     1.71956e-07    0.000340879        0.00134      
%     79       1040     1.65903e-07    0.000681758        0.00526      
%     80       1053     1.57598e-07    0.000681758        0.00525      
%     81       1066     1.57598e-07     0.00136352        0.00525      
%     82       1079     1.51097e-07    0.000340879        0.00131      
%     83       1092     1.45498e-07    0.000681758        0.00511      
%     84       1105     1.37755e-07    0.000681758        0.00509      
%     85       1118     1.30182e-07    0.000681758        0.00502      
%     86       1131     1.22817e-07    0.000681758        0.00495      
%     87       1144     1.15653e-07    0.000681758        0.00489      
%     88       1157      1.0868e-07    0.000681758        0.00482      
%     89       1170     1.01893e-07    0.000681758        0.00476      
%     90       1183     9.52847e-08    0.000681758         0.0047      
%     91       1196     8.88531e-08    0.000681758        0.00465      
%     92       1209     8.25949e-08    0.000681758        0.00459      
%     93       1222     7.65082e-08    0.000681758        0.00454      
%     94       1235      7.0592e-08    0.000681758        0.00449      
%     95       1248     6.48501e-08    0.000681758        0.00444      
%     96       1261     5.92838e-08    0.000681758        0.00439      
%     97       1274     5.38992e-08    0.000681758        0.00435      
%     98       1287     4.87015e-08    0.000681758        0.00431      
%     99       1300      4.3696e-08    0.000681758        0.00426      
%    100       1313     3.88953e-08    0.000681758        0.00424      
%    101       1326     3.43071e-08    0.000681758        0.00421      
%    102       1339     2.99433e-08    0.000681758        0.00419      
%    103       1352      2.5819e-08    0.000681758        0.00416      
%    104       1365     2.19481e-08    0.000681758        0.00414      
%    105       1378     1.83478e-08    0.000681758        0.00412      
%    106       1391     1.50344e-08    0.000681758         0.0041      
%    107       1404     1.20268e-08    0.000681758        0.00408      
%    108       1417     9.34869e-09    0.000681758        0.00406      
%    109       1430     7.01889e-09    0.000681758        0.00404      
%    110       1443     5.06243e-09    0.000681758        0.00402      
%    111       1456     4.00556e-09    0.000741534        0.00473      
%    112       1469     1.54219e-09     0.00017044       0.000256      
%    113       1482     1.09533e-09    0.000340879          0.001      
%    114       1495     1.09533e-09    0.000681758          0.001      
%    115       1508     7.77606e-10     0.00017044       0.000256      
%    116       1521     4.98278e-10    0.000340879       0.000998      
%    117       1534     2.44533e-10    0.000340879       0.000993      
%    118       1547     1.07561e-10     0.00033576       0.000957      
%    119       1560     4.06888e-12    0.000123572       0.000129      
%    120       1573     1.19259e-13    6.28433e-05       3.36e-05      
%    121       1586     6.06729e-18    4.49365e-06       1.58e-07      
% 
% Local minimum found.
% 
% Optimization completed because the size of the gradient is less than
% the value of the optimality tolerance.
% 
% <stopping criteria details>
% 
% Optimization Results (Tikhonov Regularized):
% Final Residual Norm (augmented): 2.4632e-09
% Exit Flag: 1
% 
% Parameter Recovery Accuracy (Tikhonov Regularized):
% Mean Relative Error: 0.00%
% Max Relative Error: 0.00%
% Saved figure to Parameter_Recovery_1_No_Noise,_No_Regularization.png
% Saved figure to Signal_Fit_1_No_Noise,_No_Regularization.png
% --- SCENARIO 1 COMPLETE ---