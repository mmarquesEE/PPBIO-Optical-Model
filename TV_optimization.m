% Main script to run the Tikhonov regularization scenarios
clearvars; close all; clc;

fprintf('This script will demonstrate the necessity of Tikhonov regularization for parameter estimation.\n');
fprintf('It will run three scenarios and generate plots for each.\n\n');


% --- SCENARIO 3: The Solution ---
% With noise AND Tikhonov regularization. We expect a good parameter recovery,
% showing how regularization stabilizes the solution.
fprintf('--- RUNNING SCENARIO 3: With Noise, With Regularization ---\n');
run_tv_estimation('With Noise, With TV Regularization', 0.02, 5e-3);
fprintf('--- SCENARIO 3 COMPLETE ---\n\n');

function run_tv_estimation(scenario_title, noise_level, lambda_tv)
    % Shared parameters
    gridN_x = 10; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [5,5]; ads_y_range = [2,3];
    ads_nx = (ads_x_range(2)-ads_x_range(1)+1); ads_ny = (ads_y_range(2)-ads_y_range(1)+1);
    num_ads_cells =  ads_nx * ads_ny; 
    
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

    % Generate 1 experiments
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
    true_params_vec = [true_kon(:); true_koff(:); true_smax(:)];

        % Generate "experimental" data with noise
    fprintf('Generating noisy data (Noise Level: %.0f%%, Num Experiments: %d)...\n', noise_level*100, n_exp);
    exp_data = cell(n_exp,1);
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        [~, s_obs] = run_single_experiment(gridN_x, gridN_y, gridN_z, true_kon, true_koff, true_smax, ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
        noise_std = noise_level * max(s_obs);
        exp_data{exp_idx} = s_obs + noise_std * randn(size(s_obs));
    end
    
    % --- Parameter Identification using TV Regularization ---
    fprintf('Starting parameter identification with TV Regularization (lambda=%.1e)...\n', lambda_tv);
    
    % Setup the problem structure to pass to the solver
    problem.exp_settings = exp_settings;
    problem.exp_data = exp_data;
    problem.gridN_x = gridN_x; problem.gridN_y = gridN_y; problem.gridN_z = gridN_z;
    problem.ads_x_range = ads_x_range; problem.ads_y_range = ads_y_range; problem.ads_layer = ads_layer;
    problem.D_coeff = D_coeff; problem.ru_to_m = ru_to_m; problem.n_exp = n_exp;
    problem.ads_nx = ads_nx; problem.ads_ny = ads_ny;
    problem.num_ads_cells = num_ads_cells;
    problem.homog_params = homog_params;
    problem.N_params = length(true_params_vec);

    % Solve using the new TV solver
    opt_params = solve_tv_irls(problem, lambda_tv, 10); % 10 IRLS iterations is a good start

    % --- Analyze and Plot Results ---
    fprintf('\nOptimization with TV Regularization Complete.\n');
    param_error = 100 * norm(opt_params - true_params_vec) / norm(true_params_vec);
    fprintf('Final Relative Parameter Error: %.2f%%\n', param_error);
    
    % Initial guess for plotting
    homog_kon = homog_params(1);
    homog_koff = homog_params(2);
    homog_smax_cell = homog_params(3) / num_ads_cells;
    p0 = [homog_kon * ones(num_ads_cells,1); homog_koff * ones(num_ads_cells,1); homog_smax_cell * ones(num_ads_cells,1)];
    
    plot_parameter_recovery(true_params_vec, opt_params, p0, ...
        size(true_kon), size(true_koff), size(true_smax), scenario_title);
        
    plot_signal_predictions(opt_params, exp_settings, exp_data, ...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer, ...
        D_coeff, ru_to_m, n_exp, scenario_title);
end

% ================== NEW TV SOLVER and HELPERS ==================
function opt_params = solve_tv_irls(problem, lambda, max_irls_iters)
    % Solves the TV-regularized problem using Iteratively Reweighted Least Squares.
    p = problem;
    
    % --- Get Initial Guess and Bounds ---
    homog_kon = p.homog_params(1); 
    homog_koff = p.homog_params(2);
    homog_smax_cell = p.homog_params(3) / p.num_ads_cells;
    p0 = [log10(homog_kon * ones(p.num_ads_cells,1)); ...
          log10(homog_koff * ones(p.num_ads_cells,1)); ...
          log10(homog_smax_cell * ones(p.num_ads_cells,1))];
    lb = [repmat(log10(1e2), p.num_ads_cells, 1); repmat(log10(1e-4), p.num_ads_cells, 1); repmat(log10(1e-3), p.num_ads_cells, 1)];
    ub = [repmat(log10(1e5), p.num_ads_cells, 1); repmat(log10(1e0), p.num_ads_cells, 1); repmat(log10(1e1), p.num_ads_cells, 1)];
    
    % --- Setup for IRLS ---
    L = create_spatial_regularization_operator(p.ads_nx, p.ads_ny);
    optim_opts = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'Display', 'iter', 'UseParallel', true, 'MaxIterations', 100, 'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8);
    
    log_p_k = p0; % Start with the initial guess
    epsilon = 1e-6; % Small stabilization parameter

    for k = 1:max_irls_iters
        fprintf('  IRLS Iteration %d/%d\n', k, max_irls_iters);
        
        % 1. Calculate weights based on the current parameter estimate
        Lp = L * log_p_k;
        weights = 1 ./ sqrt(Lp.^2 + epsilon);
        W = spdiags(weights, 0, length(weights), length(weights));
        
        % 2. Define the *weighted* L2 residual function
        residual_fun_weighted = @(log_params) ...
            compute_weighted_residuals(log_params, p, lambda, sqrt(W)*L);
            
        % 3. Solve the weighted least-squares problem for the next iterate
        log_p_k = lsqnonlin(residual_fun_weighted, log_p_k, lb, ub, optim_opts);
    end
    
    opt_params = 10.^log_p_k;
end

function residuals_aug = compute_weighted_residuals(log_params, problem, lambda, L_weighted)
    % This function computes the residual for a weighted L2 problem,
    % which is a single step inside the IRLS algorithm for TV.
    p = problem;
    
    params = 10.^log_params;
    kon_ads = reshape(params(1:p.num_ads_cells), [p.ads_nx, p.ads_ny]);
    koff_ads = reshape(params(p.num_ads_cells+1:2*p.num_ads_cells), [p.ads_nx, p.ads_ny]);
    smax_ads = reshape(params(2*p.num_ads_cells+1:end), [p.ads_nx, p.ads_ny]);
    
    model_residuals = [];
    sim_data = cell(p.n_exp, 1);

    % This loop can be parallelized with 'parfor' if you have the toolbox
    for i = 1:p.n_exp
        setting = p.exp_settings(i);
        [~, s_sim] = run_single_experiment(...
            p.gridN_x, p.gridN_y, p.gridN_z, kon_ads, koff_ads, smax_ads,...
            p.ads_x_range, p.ads_y_range, p.ads_layer, setting, p.D_coeff, p.ru_to_m);
        sim_data{i} = s_sim;
    end
    
    for i = 1:p.n_exp
        model_residuals = [model_residuals; (sim_data{i} - p.exp_data{i})];
    end

    % Add the weighted regularization term
    reg_term_vector = sqrt(lambda) * (L_weighted * log_params(:));
    residuals_aug = [model_residuals; reg_term_vector];
end

function L_full = create_spatial_regularization_operator(nx, ny)
    % Creates a sparse finite difference operator for spatial regularization.
    N = nx * ny;
    e_x = ones(N, 1);
    Dx = spdiags([-e_x, e_x], [0, 1], N, N);
    for i = 1:(ny-1)
        Dx(i*nx, i*nx + 1) = 0;
    end
    Dx(ny*nx, :) = [];
    e_y = ones(N, 1);
    Dy = spdiags([-e_y, e_y], [0, nx], N, N);
    Dy(N-nx+1:end, :) = [];
    L_space = [Dx; Dy];
    L_full = blkdiag(L_space, L_space, L_space);
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
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    
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