% Main script to run the Tikhonov regularization scenarios
clearvars; close all; clc;

% main_run_lcurve_analysis.m
% This script first finds the optimal lambda via L-curve analysis and then
% runs the final, regularized parameter estimation.

clearvars; close all; clc;

fprintf('--- L-CURVE REGULARIZATION ANALYSIS ---\n\n');

% =========================================================================
% STEP 1: SETUP AND DATA GENERATION (Same as before)
% =========================================================================
fprintf('STEP 1: Setting up the problem and generating noisy data...\n');
% Shared parameters
gridN_x = 10; gridN_y = 5; gridN_z = 3;
ads_layer = 1;
ads_x_range = [5,5]; ads_y_range = [2,3];

ads_nx = ads_x_range(2) - ads_x_range(1) + 1;
ads_ny = ads_y_range(2) - ads_y_range(1) + 1;
num_ads_cells = ads_nx * ads_ny;

% Homogeneous parameters for initial guess
homog_params = [9.4e3, 0.0078, 1.0]; % kon, koff, smax_total

% Fixed parameters
D_coeff = 6e-3;
ru_to_m = 1e-6;
base_max_velocity = 8.3;
noise_level = 0.02; % 5% noise

% Experiment settings
n_exp = 3;
T1 = 2000; T2 = 4000; T3 = 3*T1;
t_total = 4*T1;
c_diss = 0; c1 = 3.3e-6; c2 = 1e-6;
exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);

% Create and extract true heterogeneous parameters
[kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
    create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
true_kon = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
true_koff = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
true_smax = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
true_params_vec = [true_kon(:); true_koff(:); true_smax(:)];
N_params = length(true_params_vec);

% Generate "experimental" data with noise
exp_data = cell(n_exp,1);
parfor exp_idx = 1:n_exp
    setting = exp_settings(exp_idx);
    [~, s_obs] = run_single_experiment(...
        gridN_x, gridN_y, gridN_z, true_kon, true_koff, true_smax,...
        ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m);
    noise_std = noise_level * max(s_obs);
    exp_data{exp_idx} = s_obs + noise_std * randn(size(s_obs));
end
fprintf('Data generation complete.\n\n');

% =========================================================================
% STEP 2: L-CURVE ANALYSIS TO FIND OPTIMAL LAMBDA
% =========================================================================
fprintf('STEP 2: Performing L-curve analysis to find optimal lambda...\n');

% Define range of lambda values to test
lambda_values = logspace(-4, 2, 5);

% Define problem structure to pass to the function
problem.exp_settings = exp_settings;
problem.exp_data = exp_data;
problem.gridN_x = gridN_x; problem.gridN_y = gridN_y; problem.gridN_z = gridN_z;
problem.ads_x_range = ads_x_range; problem.ads_y_range = ads_y_range; problem.ads_layer = ads_layer;
problem.D_coeff = D_coeff; problem.ru_to_m = ru_to_m;
problem.n_exp = n_exp;
problem.ads_nx = ads_nx; problem.ads_ny = ads_ny;
problem.num_ads_cells = num_ads_cells;
problem.homog_params = homog_params;
problem.N_params = N_params;

% Run the L-curve analysis
optimal_lambda = find_optimal_lambda_lcurve(problem, lambda_values);

fprintf('L-curve analysis complete. Optimal lambda found: %.4e\n\n', optimal_lambda);

% =========================================================================
% STEP 3: RUN FINAL ESTIMATION WITH OPTIMAL LAMBDA
% =========================================================================
fprintf('STEP 3: Running final estimation with the optimal lambda.\n');

scenario_title = sprintf('Optimal Regularization (lambda=%.2e)', optimal_lambda);
run_estimation_scenario(scenario_title, noise_level, optimal_lambda);

fprintf('--- ANALYSIS COMPLETE ---\n');

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

% find_optimal_lambda_lcurve.m
function optimal_lambda = find_optimal_lambda_lcurve(problem, lambda_values)
    % Unpack problem structure
    p = problem; 
    
    % Get initial guess and bounds
    homog_kon = p.homog_params(1);
    homog_koff = p.homog_params(2);
    homog_smax_cell = p.homog_params(3) / p.num_ads_cells;
    p0 = [log10(homog_kon * ones(p.num_ads_cells,1)); ...
          log10(homog_koff * ones(p.num_ads_cells,1)); ...
          log10(homog_smax_cell * ones(p.num_ads_cells,1))];
    lb = [repmat(log10(1e2), p.num_ads_cells, 1); repmat(log10(1e-4), p.num_ads_cells, 1); repmat(log10(1e-3), p.num_ads_cells, 1)];
    ub = [repmat(log10(1e5), p.num_ads_cells, 1); repmat(log10(1e0), p.num_ads_cells, 1); repmat(log10(1e1), p.num_ads_cells, 1)];
    
    % Create the spatial regularization operator
    L_matrix = create_spatial_regularization_operator(p.ads_nx, p.ads_ny);
    log_params_ref = zeros(p.N_params, 1);

    % Optimization options
    optim_opts = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', 'UseParallel', true, 'FiniteDifferenceType', 'central',...
        'MaxIterations', 15, 'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8);

    % Preallocate results
    rho = zeros(length(lambda_values), 1); % Residual norm ||Ax-b||^2
    eta = zeros(length(lambda_values), 1); % Solution norm ||Lx||^2
    
    fprintf('Testing %d lambda values...\n', length(lambda_values));
    tic;
    
    % Loop over all lambda values
    for i = 1:length(lambda_values)
        lambda = lambda_values(i);
        
        % Define the residual function for this lambda
        residual_fun = @(log_params) compute_residuals_log_tikhonov(...
            log_params, p.exp_settings, p.exp_data, ...
            p.gridN_x, p.gridN_y, p.gridN_z, ...
            p.ads_x_range, p.ads_y_range, p.ads_layer, p.D_coeff, p.ru_to_m, ...
            p.n_exp, lambda, log_params_ref, L_matrix);
        
        % Run optimization
        [opt_log_params, resnorm] = lsqnonlin(residual_fun, p0, lb, ub, optim_opts);
        
        % Deconstruct the augmented residual norm
        num_model_residuals = sum(cellfun(@numel, p.exp_data));
        num_reg_residuals = size(L_matrix, 1);
        
        % The total resnorm = rho + eta
        % We need to separate them.
        residuals_unaugmented = residual_fun(opt_log_params);
        residuals_model = residuals_unaugmented(1:num_model_residuals);
        residuals_reg = residuals_unaugmented(num_model_residuals+1:end);

        rho(i) = norm(residuals_model)^2;
        % eta = ||L*p||^2 = (||sqrt(lambda)*L*p||^2) / lambda
        eta(i) = norm(residuals_reg)^2 / lambda;
        
        fprintf('  Lambda = %.2e, rho = %.2e, eta = %.2e\n', lambda, rho(i), eta(i));
    end
    
    toc;
    
    % Find the corner of the L-curve
    [k, idx] = max(calculate_curvature(log10(rho), log10(eta)));
    optimal_lambda = lambda_values(idx);
    
    % Plot the L-curve
    figure('Name', 'L-Curve Analysis', 'Position', [200, 200, 700, 600]);
    loglog(rho, eta, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    hold on;
    loglog(rho(idx), eta(idx), 'r-s', 'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');
    xlabel('Residual Norm ||\Phi(p) - s_{obs}||_2^2');
    ylabel('Solution Norm ||L p||_2^2');
    title('L-Curve: Trade-off between Data Fit and Solution Smoothness');
    legend('L-Curve', sprintf('Optimal \lambda = %.2e', optimal_lambda), 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);
    
    fig_filename = 'L_Curve_Analysis.png';
    saveas(gcf, fig_filename);
    fprintf('Saved L-curve plot to %s\n', fig_filename);
end

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

% calculate_curvature.m
function [curvature, max_idx] = calculate_curvature(x, y)
    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Use finite differences to approximate derivatives
    dx = gradient(x);
    dy = gradient(y);
    ddx = gradient(dx);
    ddy = gradient(dy);
    
    % Curvature formula for a parametric curve
    num = dx .* ddy - dy .* ddx;
    den = (dx.^2 + dy.^2).^(3/2);
    
    curvature = abs(num ./ den);
    
    % The best point is usually not at the very ends, so we might
    % ignore the first and last points.
    [~, max_idx] = max(curvature(2:end-1));
    max_idx = max_idx + 1; % Adjust index
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
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
    
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
        num_points = 500;
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