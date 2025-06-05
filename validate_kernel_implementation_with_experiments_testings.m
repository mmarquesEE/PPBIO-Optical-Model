function validate_kernel_implementation_with_experiments()
    clearvars; close all; clc;
    % Shared parameters
    gridN_x = 10; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [5,6]; ads_y_range = [2,2];
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Homogeneous parameters
    homog_params = [9.4e3, 0.0078, 1.0]; % kon, koff, smax_total
    smax_per_cell = homog_params(3) / num_ads_cells;
    
    % Fixed parameters
    D_coeff = 6e-3;
    ru_to_m = 1e-6;
    base_max_velocity = 8.3;
    
    % Example pulse parameters
    
    % In your main validation function:
    base_max_velocity = 8.3;  % µm/s
    T_total = 3200;            % Total experiment duration (s)
    c1 = 3.3e-6;               % Primary concentration (M)
    c2 = 1e-6;                 % Secondary concentration (M)
    
    % Calculate number of experiments (M >= 3⌈log₃n⌉)
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    n_params = 3 * num_ads_cells;
    M = max(3, 3 * ceil(log(n_params)/log(3)));
    
    % Generate experiments
    exp_settings = generate_experiments(M, base_max_velocity, T_total, c1, c2);
    n_exp = length(exp_settings);
    % Display generated experiments
    parfor i = 1:length(exp_settings)
        fprintf('Experiment %d:\n', i);
        fprintf('  Concentrations: [%.2e, %.2e, %.2e]\n', exp_settings(i).pulse_concs);
        fprintf('  Velocity: %.2f µm/s\n', exp_settings(i).max_velocity);
    end

    
    % Create heterogeneous parameters
    [kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
        create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
    % Preallocate parfor identifiability analysis
    s_obs_final_base = zeros(n_exp,1); % Final s_obs values parfor baseline (each experiment)
    
    parfor exp_idx = 1:n_exp
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
        
        % Initial condition parfor s (empty)
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
        
        % Extract adsorption region parfor state variables
        K_homog_ads = K_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        Q_homog_ads = Q_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        s_homog_ads = s_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        
        K_heterog_ads = K_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        Q_heterog_ads = Q_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        s_heterog_ads = s_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
        
        % Get initial surface concentration (t=0)
        s0_homog_ads = squeeze(s_homog_ads(1, :, :, :));
        s0_heterog_ads = squeeze(s_heterog_ads(1, :, :, :));
        
        % Compute alpha (kon * smax) parfor each cell
        alpha_homog_ads = kon_homog_ads .* smax_homog_ads;
        alpha_heterog_ads = kon_heterog_ads .* smax_heterog_ads;
        
        % Reshape parfor broadcasting
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
        
        % Direct s_obs parfor verification
        s_obs_homog_direct = compute_s_obs(s_homog, ads_x_range, ads_y_range, ads_layer);
        s_obs_heterog_direct = compute_s_obs(s_heterog, ads_x_range, ads_y_range, ads_layer);
        
        % Store final s_obs value parfor identifiability analysis
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
        fprintf('Computing Jacobian parfor Experiment %d...\n', exp_idx);
        
        % Get time vector and composite signal parfor this experiment
        [~, s_obs] = run_single_experiment(...
            gridN_x, gridN_y, gridN_z, kon_ads, koff_ads, smax_ads,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        
        % Compute Jacobian via finite differences
        J_exp = zeros(length(s_obs), N_params);
        h = 1e-5;  % Finite difference step
        
        for k = 1:N_params
            p_pert = p_true;
            p_pert(k) = p_pert(k) + h;
            
            % Split parameters into kon/koff/smax components
            kon_pert = reshape(p_pert(1:numel(kon_ads)), size(kon_ads));
            koff_pert = reshape(p_pert(numel(kon_ads)+1:numel(kon_ads)+numel(koff_ads)), size(koff_ads));
            smax_pert = reshape(p_pert(numel(kon_ads)+numel(koff_ads)+1:end), size(smax_ads));
            
            % Run simulation with perturbed parameters
            [~, s_pert] = run_single_experiment(...
                gridN_x, gridN_y, gridN_z, kon_pert, koff_pert, smax_pert,...
                ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
            
            % Compute derivative
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
    
    % Add text annotations
    annotation('textbox', [0.15, 0.15, 0.3, 0.1], 'String', ...
        sprintf('Rank = %d/%d\nCond = %.1e', rankJ, N_params, cond_number), ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white');
    
    % ============== ENHANCED PARAMETER IDENTIFICATION ================
    fprintf('\nStarting Enhanced Parameter Identification...\n');
    
    % Extract true heterogeneous parameters
    true_kon = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    true_koff = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    true_smax = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    % Generate "experimental" data from true parameters
    exp_data = cell(n_exp,1);
    parfor exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [t, s_obs] = run_single_experiment(...
            gridN_x, gridN_y, gridN_z, true_kon, true_koff, true_smax,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        noise_level = 0.02;  % 1% noise
        exp_data{exp_idx} = s_obs.* (1 + noise_level*randn(size(s_obs)));
    end
    
    % Prepare logarithmic parameterization
    num_ads_cells = numel(true_kon);
    homog_kon = homog_params(1);
    homog_koff = homog_params(2);
    homog_smax_cell = homog_params(3) / num_ads_cells;
    
    % Initial guess in log space
    rng(2); % parfor reproducible initialization
    p0_kon = log10(homog_kon) + 5*randn(num_ads_cells,1);
    p0_koff = log10(homog_koff) + 5*randn(num_ads_cells,1);
    p0_smax = log10(homog_smax_cell) + 5*randn(num_ads_cells,1);
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
    
    % Optimization options
    optim_opts = optimoptions('lsqnonlin', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', ...
        'UseParallel',true,...
        'MaxIterations', 15, ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'FiniteDifferenceType', 'central');
    
    % Residual function
    residual_fun = @(log_params) compute_residuals_log(...
        log_params, exp_settings, exp_data, gridN_x, gridN_y, gridN_z, ...
        ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,n_exp);
    
    % Run optimization in log space
    [opt_log_params, ~, residual, exitflag] = lsqnonlin(...
        residual_fun, p0, lb, ub, optim_opts);
    
    % Convert optimized parameters back to linear space
    opt_params = 10.^opt_log_params;
    
    % Split into parameter groups
    opt_kon = opt_params(1:num_ads_cells);
    opt_koff = opt_params(num_ads_cells+1:2*num_ads_cells);
    opt_smax = opt_params(2*num_ads_cells+1:end);
    
    % Reshape to match adsorption region
    opt_kon = reshape(opt_kon, size(true_kon));
    opt_koff = reshape(opt_koff, size(true_koff));
    opt_smax = reshape(opt_smax, size(true_smax));
    
    % True parameters as vector parfor comparison
    true_params_vec = [true_kon(:); true_koff(:); true_smax(:)];
    init_params = 10.^p0; % Initial guess in linear space
    
    % Analyze results
    fprintf('\nOptimization Results:\n');
    fprintf('Final Residual Norm: %.4e\n', norm(residual));
    fprintf('Exit Flag: %d\n', exitflag);
    
    % Calculate parameter errors
    param_errors = abs(opt_params - true_params_vec) ./ true_params_vec;
    fprintf('\nParameter Recovery Accuracy:\n');
    fprintf('Mean Relative Error: %.2f%%\n', 100*mean(param_errors));
    fprintf('Max Relative Error: %.2f%%\n', 100*max(param_errors));
    
    % Plot parameter recovery
    plot_parameter_recovery(true_params_vec, opt_params, init_params, ...
        size(true_kon), size(true_koff), size(true_smax));
    
    % Plot predicted vs "experimental" signals
    plot_signal_predictions([opt_kon(:); opt_koff(:); opt_smax(:)], ...
        exp_settings, exp_data, gridN_x, gridN_y, gridN_z, ...
        ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,n_exp);

end

function exp_settings = generate_experiments(M, base_vel, T_total, c1, c2)
    % Inputs:
    %   M - Number of experiments
    %   base_vel - Base velocity (µm/s)
    %   T_total - Total experiment duration
    %   c1, c2 - Base concentrations
    
    % Non-negative concentration design via NMF
    C_pos = design_nonnegative_concentrations(M, c1, c2);
    
    % Velocity range (exponential spacing)
    v_min = base_vel * exp(-sqrt(M));
    v_max = base_vel * exp(sqrt(M));
    velocities = exp(linspace(log(v_min), log(v_max), M));
    
    % Common pulse times (equally spaced)
    pulse_times = [T_total/4, T_total/2, 3*T_total/4];
    
    % Build experiments
    exp_settings = struct();
    parfor m = 1:M
        exp_settings(m).pulse_times = pulse_times;
        exp_settings(m).pulse_concs = C_pos(m, :);
        exp_settings(m).t_total = T_total;
        exp_settings(m).max_velocity = velocities(m);
        exp_settings(m).c_diss = 0; % Dissociation buffer
    end
end


% ================== NEW HELPER FUNCTIONS ==================
function residuals = compute_residuals_log(log_params, exp_settings, exp_data, ...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,n_exp)

    % Convert from log to linear scale
    params = 10.^log_params;
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Split parameters
    kon_ads = reshape(params(1:num_ads_cells), [ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1]);
    koff_ads = reshape(params(num_ads_cells+1:2*num_ads_cells), size(kon_ads));
    smax_ads = reshape(params(2*num_ads_cells+1:end), size(kon_ads));
    
    residuals = [];
    parfor exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [~, s_sim] = run_single_experiment(...
            nx, ny, nz, kon_ads, koff_ads, smax_ads,...
            ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        
        % Stack residuals from all experiments
        residuals = [residuals; (s_sim - exp_data{exp_idx})];
    end
end
function C_pos = design_nonnegative_concentrations(M, c1, c2)
    % Create non-negative orthogonal basis
    angles = linspace(0, pi/2, max(M, 3)); % Ensure coverage
    C = zeros(M, 3);
    
    for m = 1:M
        % Use trigonometric functions parfor non-negative values
        theta = angles(mod(m-1, length(angles)) + 1);
        phi = angles(ceil(m/2));
        
        % Pulse 1: Main binding pulse
        C(m,1) = c1 * (0.2 + 0.8*cos(theta)^2);
        
        % Pulse 2: Optional secondary binding
        if mod(m,3) == 0
            C(m,2) = c2 * (0.3 + 0.7*sin(phi)^2);
        else
            C(m,2) = 0; % Dissociation phase
        end
        
        % Pulse 3: Competitive binding
        C(m,3) = c2 * (0.2 + 0.8*sin(theta)^2);
    end
    
    % Ensure minimum concentration
    min_conc = 1e-9; % 1 nM minimum
    C_pos = max(C, min_conc);
    
    % Orthogonalize while preserving non-negativity
    [U,~,~] = svd(C_pos);
    C_pos = abs(U(:,1:3)) * diag([c1, c2, c2]);
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

function plot_parameter_recovery(true_params, opt_params, initial_guess, sz_kon, sz_koff, sz_smax)
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
    
    sgtitle('Parameter Recovery Results');
end

function plot_signal_predictions(opt_params, exp_settings, exp_data, ...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,n_exp)
    
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    kon_ads = reshape(opt_params(1:num_ads_cells), [ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1]);
    koff_ads = reshape(opt_params(num_ads_cells+1:2*num_ads_cells), size(kon_ads));
    smax_ads = reshape(opt_params(2*num_ads_cells+1:end), size(kon_ads));
    
    figure('Position', [100, 100, 1200, 800]);
    parfor exp_idx = 1:n_exp
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
    sgtitle('Signal Prediction vs Experimental Data');
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
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    
    % Preallocate results
    t_all = [];
    y_all = [];
    
    % Process each time segment
    num_segments = length(t_breaks) - 1;
    for seg = 1:num_segments
        t_start = t_breaks(seg);
        t_end = t_breaks(seg+1);
        c0_seg = concentrations(seg);
        
        % Determine time points parfor segment
        num_points = max(10, ceil(200 * (t_end - t_start) / (t_breaks(end) - t_breaks(1))));
        tspan = linspace(t_start, t_end, num_points);
        
        % Run simulation parfor segment
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
        
        % Update initial condition parfor next segment
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