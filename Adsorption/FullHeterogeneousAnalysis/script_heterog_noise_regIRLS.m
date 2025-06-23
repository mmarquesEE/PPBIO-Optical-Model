%% Full Heterogeneous Model (with Line Plot Visualization)
clearvars; close all; clc;

%% --- STEP 0.1: Parameters Definition ---
% ========================================================================
generate_video_frames = false; 
gridN_x = 22; gridN_y = 5; gridN_z = 3; ads_layer = 1;
ads_x_range = [10,11]; ads_y_range = [1,3];

% --- Key Dimensions ---
ads_x_dim = ads_x_range(2) - ads_x_range(1) + 1;
ads_y_dim = ads_y_range(2) - ads_y_range(1) + 1;
num_sites = ads_x_dim * ads_y_dim; % Total number of heterogeneous sites

% Homogeneous parameters & Fixed parameters
homog_params = [9.4e3, 0.0078, 2960];
D_coeff = 6e-5; ru_to_m = 1e-10; base_max_velocity = 8.3;
grid_size_x = 11.0; grid_size_z = 0.3;
dx = grid_size_x / gridN_x; dz = grid_size_z / gridN_z;
n_exp = 1; T1 = 2200; T2 = 2*T1; T3 = 3*T1; t_total = 4*T1;
c_diss = 0; c1 = 3.3e-4; c2 = 0.21e-4;
exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);
model_config.gridN_x = gridN_x; model_config.gridN_y = gridN_y; model_config.gridN_z = gridN_z;
model_config.dx = dx; model_config.dz = dz; model_config.ads_x_range = ads_x_range;
model_config.ads_y_range = ads_y_range; model_config.ads_layer = ads_layer;
model_config.D_coeff = D_coeff; model_config.ru_to_m = ru_to_m;

% Create ground truth using the full heterogeneity function
[kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
    create_ground_truth_heterogeneity_full(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);

%=========================================================================
%% --- STEP 1: GENERATE "EXPERIMENTAL" DATA ---
%=========================================================================
fprintf('\nGenerating ground-truth data...\n');
tic;
exp_data = cell(n_exp, 1);
total_rows = 0;
for exp_idx = 1:n_exp
    setting = exp_settings(exp_idx);
    velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    s0_grid = zeros(gridN_x, gridN_y, gridN_z);
    
    [t_exp, ~, s_heterog_exp] = simulate_3d_flow_model_with_pulses(gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid,dx,dz);
    
    s_obs_by_line_clean = compute_s_obs_by_line(s_heterog_exp, ads_x_range, ads_y_range, ads_layer);
    
    noise_level = 0.01;
    max_response = max(s_obs_by_line_clean(:));
    noise_std_dev = noise_level * max_response;
    additive_noise = noise_std_dev * randn(size(s_obs_by_line_clean));
    s_obs_by_line_noisy = s_obs_by_line_clean + additive_noise;
    
    data_struct.time = t_exp;
    data_struct.signals_clean = s_obs_by_line_clean;
    data_struct.signals = s_obs_by_line_noisy;
    exp_data{exp_idx} = data_struct;
    total_rows = total_rows + numel(exp_data{exp_idx}.signals_clean);
end
fprintf('Data generation complete.\n'); toc;

%=========================================================================
%% --- STEP 2: IDENTIFIABILITY & PARAMETER IDENTIFICATION (REGULARIZED) ---
%=========================================================================
fprintf('\nStarting 2D Regularized Parameter Identification...\n');

% Create the "true" 2D parameter vector (no change here)
kon_ads_2D = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
koff_ads_2D = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
smax_ads_2D = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
p_true_2D = [kon_ads_2D(:); koff_ads_2D(:); smax_ads_2D(:)];
N_params_2D = length(p_true_2D);

% Setup initial guess and bounds (no change here)
homog_kon = homog_params(1);
homog_koff = homog_params(2);
smax_per_site = homog_params(3) / num_sites;
p0_base = [log10(homog_kon) * ones(num_sites, 1); ...
           log10(homog_koff) * ones(num_sites, 1); ...
           log10(smax_per_site) * ones(num_sites, 1)];
lb = [repmat(log10(1e3), num_sites, 1); repmat(log10(1e-6), num_sites, 1); repmat(log10(1e-4), num_sites, 1)];
ub = [repmat(log10(1e7), num_sites, 1); repmat(log10(1e-1), num_sites, 1); repmat(log10(smax_per_site * 10), num_sites, 1)];
rng('default');
p0_rand = p0_base + 0.1 * randn(size(p0_base));
p0 = max(lb, min(ub, p0_rand));
p_current_log = p0; % This is our initial guess


% ============================================================
% --- IRLS Setup ---
% ============================================================
% 1. Regularization and IRLS parameters
lambda = 0.2;          % The overall regularization strength (tune this)
irls_iterations = 10;   % Number of outer loops (5-15 is typical)
epsilon = 1e-6;         % Small constant to avoid division by zero

% 2. Create the regularization matrix L (this is done only once)
L = create_regularization_matrix(ads_x_dim, ads_y_dim);

% 3. Initialize the weight matrix W as the identity matrix
W = speye(size(L, 1));

% 4. Set options for the INNER lsqnonlin calls (fewer iterations needed)
optim_plot_fun = @(log_p, optim_v, state) optimPlotter_2D_as_lines(log_p, optim_v, state, p_true_2D, exp_data, model_config);

% Set optimization options (no change here)
optim_opts_inner  = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', 'Display', 'iter', 'MaxIterations', 10, ...
    'UseParallel', true, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OutputFcn', optim_plot_fun);
    
% ============================================================
% --- The Main IRLS Loop ---
% ============================================================
tic;
for iter = 1:irls_iterations
    fprintf('\n--- Starting IRLS Iteration %d of %d ---\n', iter, irls_iterations);
    
    % Define the residual function for this iteration, using the current W
    residual_fun = @(log_params) compute_weighted_regularized_residuals(log_params, W, L, lambda, exp_settings, exp_data, model_config);
    
    % --- Run the inner optimization ---
    % We use p_current_log (the result from the last iteration) as the new starting guess
    [opt_log_params, resnorm, residual, exitflag, output, lambda_lsq, J_opt] = ...
        lsqnonlin(residual_fun, p_current_log, lb, ub, optim_opts_inner);
    
    % --- Update the solution for the next iteration ---
    p_current_log = opt_log_params;
    p_current_linear = 10.^p_current_log;
    
    % --- Update the weights for the next iteration ---
    Lp = L * p_current_linear;
    new_weights_vector = 1 ./ (sqrt(Lp.^2 + epsilon)); % Using sqrt(x^2) is more robust than abs(x)
    W = spdiags(new_weights_vector, 0, size(L,1), size(L,1));
    
    fprintf('Finished IRLS Iteration %d. Norm of residual: %.4e\n', iter, resnorm);
end
toc;

% --- Analyze and plot final results ---
p_init = 10.^p0;
plot_final_2D_results_full(p_true_2D, opt_log_params, p_init, resnorm, residual, J_opt, output, exp_data, model_config);

