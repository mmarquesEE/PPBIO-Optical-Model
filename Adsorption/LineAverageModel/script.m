%% Line Average Model
clearvars; close all; clc;

% --- Create the EPS folder if it doesn't exist ---
eps_folder_path = 'Adsorption/LineAverageModel/Figures/EPS';
if ~exist(eps_folder_path, 'dir')
    mkdir(eps_folder_path);
end

%% --- STEP 0.1: Parameters Definition ---
% ========================================================================
generate_video_frames = true; % Mude para 'false' para pular a criação do vídeo e acelerar o script
gridN_x = 22; gridN_y = 5; gridN_z = 3;ads_layer = 1;
ads_x_range = [5,15]; ads_y_range = [1,5];
% Get number of lines in adsorption region
ads_y_dim = ads_y_range(2) - ads_y_range(1) + 1;
% Homogeneous parameters
homog_params = [9.4e3, 0.0078, 2960]; % kon, koff, smax_total
% Fixed parameters
D_coeff = 6e-5;ru_to_m = 1e-10;base_max_velocity = 8.3;
grid_size_x = 11.0; % mm
grid_size_z = 0.3; % mm
dx = grid_size_x / gridN_x; % Size of one grid cell in x-direction (mm)
dz = grid_size_z / gridN_z; % Size of one grid cell in z-direction (mm)
% Example pulse parameters
n_exp = 1;T1 = 2200; T2 = 2*T1; T3 = 3*T1;t_total = 4*T1;
c_diss = 0;c1 = 3.3e-4;c2 = 0.21e-4;
exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);
model_config.gridN_x = gridN_x;
model_config.gridN_y = gridN_y;
model_config.gridN_z = gridN_z;
model_config.dx = dx; model_config.dz = dz;
model_config.ads_x_range = ads_x_range;
model_config.ads_y_range = ads_y_range;
model_config.ads_layer = ads_layer;
model_config.D_coeff = D_coeff;
model_config.ru_to_m = ru_to_m;
% Create ground truth 2D heterogeneous parameters
[kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
    create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
%=========================================================================
%% --- STEP 0.2: Visualizing the surfaces ---
% ========================================================================
fprintf('\nVisualizing the surfaces as 2D heatmaps...\n');
% --- Define common plotting properties for publication ---
target_fig_width_cm = 8.4; % IEEE single column width
base_font_size = 8;        % Match your paper's caption font size (e.g., 8pt)
line_width = 1.0;          % Thinner lines for smaller figures
% Extract the 2D slice of each parameter for the entire grid
kon_slice = squeeze(kon_grid_heterog(:, :, ads_layer));
koff_slice = squeeze(koff_grid_heterog(:, :, ads_layer));
smax_slice = squeeze(smax_grid_heterog(:, :, ads_layer));
% Create the figure and subplots
fig1 = figure('Name', 'Surface Parameter Heatmaps', 'Units', 'centimeters');
% --- Define the rectangle's position and size from your range variables ---
% Position format: [x_start, y_start, width, height]
% We subtract 0.5 to center the rectangle around the pixels.
rect_pos = [ads_x_range(1)-0.5, ads_y_range(1)-0.5, ...
            ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1];
            
% --- Define the explicit axes for the plot ---
x_axis_coords = 1:gridN_x;
y_axis_coords = 1:gridN_y;

% --- Create a tiled layout with minimal spacing ---
% Use 'compact' or 'tight' to remove unnecessary padding and space between tiles.
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- 2D Plot for k_on ---
ax1 = nexttile;
imagesc(ax1, x_axis_coords, y_axis_coords, kon_slice');
axis(ax1, 'xy'); 
hold(ax1, 'on');
rectangle(ax1, 'Position', rect_pos, 'EdgeColor', 'r', 'LineWidth', line_width, 'LineStyle', '--');
hold(ax1, 'off');
title(ax1, 'k_{on}', 'FontSize', base_font_size);
ylabel(ax1, 'Line Index (y)', 'FontSize', base_font_size);
colorbar(ax1);

% --- 2D Plot for k_off ---
ax2 = nexttile;
imagesc(ax2, x_axis_coords, y_axis_coords, koff_slice');
axis(ax2, 'xy');
hold(ax2, 'on');
rectangle(ax2, 'Position', rect_pos, 'EdgeColor', 'r', 'LineWidth', line_width, 'LineStyle', '--');
hold(ax2, 'off');
title(ax2, 'k_{off}', 'FontSize', base_font_size);
xlabel(t, 'Position along Flow (x)', 'FontSize', base_font_size); % Add a shared x-label
set(ax2, 'YTickLabel', []); % Remove redundant Y-axis labels

colorbar(ax2);

% --- 2D Plot for s_max ---
ax3 = nexttile;
imagesc(ax3, x_axis_coords, y_axis_coords, smax_slice');
axis(ax3, 'xy');
hold(ax3, 'on');
rectangle(ax3, 'Position', rect_pos, 'EdgeColor', 'r', 'LineWidth', line_width, 'LineStyle', '--');
hold(ax3, 'off');
title(ax3, 's_{max}', 'FontSize', base_font_size);
set(ax3, 'YTickLabel', []); % Remove redundant Y-axis labels

colorbar(ax3);

% --- Set font size for all axes in the layout ---
set([ax1, ax2, ax3], 'FontSize', base_font_size - 1);

% --- SAVE THE ENTIRE FIGURE ---
% Use the helper function to set the final size to 8.4cm and save
save_pub_fig(fig1, 'Adsorption/LineAverageModel/Figures/figure_1_surf_params_2D_combined', target_fig_width_cm);
close(fig1); % Close figure after saving

% --- SAVE EACH SUBPLOT INDIVIDUALLY ---
% fprintf('Saving individual surface parameter heatmaps as EPS files...\n');
% base_path = 'Adsorption/LineAverageModel/Figures/EPS/';
% save_subplot_as_eps(ax1, [base_path, 'surf_params_heatmap_kon.eps']);
% save_subplot_as_eps(ax2, [base_path, 'surf_params_heatmap_koff.eps']);
% save_subplot_as_eps(ax3, [base_path, 'surf_params_heatmap_smax.eps']);
% fprintf('Finished saving individual surface parameter heatmaps.\n');
% =========================================================================
%% --- STEP 1: GENERATE "EXPERIMENTAL" DATA ---
% By creating this data upfront, we know the exact size of all outputs
% and can reuse the clean data later.
% =========================================================================
fprintf('\nGenerating ground-truth data for all experiments...\n');
tic;
exp_data = cell(n_exp, 1);
total_rows = 0;
for exp_idx = 1:n_exp
    setting = exp_settings(exp_idx);
    velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    s0_grid = zeros(gridN_x, gridN_y, gridN_z);
    
    [t_exp, ~, s_heterog_exp] = simulate_3d_flow_model_with_pulses(...
        gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
        velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid,dx,dz);
    
    s_obs_by_line_clean = compute_s_obs_by_line(s_heterog_exp, ads_x_range, ads_y_range, ads_layer);
    
    noise_level = 0.01; % 1% de ruído, como antes
        
    % 1. Calcule a resposta máxima no sinal limpo para definir a escala do ruído
    max_response = max(s_obs_by_line_clean(:));
    
    % 2. Defina o desvio padrão do ruído como 1% da resposta máxima
    noise_std_dev = noise_level * max_response;
    
    % 3. Gere o ruído aditivo com desvio padrão constante
    additive_noise = noise_std_dev * randn(size(s_obs_by_line_clean));
    
    % 4. Adicione o ruído ao sinal limpo
    s_obs_by_line_noisy = s_obs_by_line_clean + additive_noise;
    
    % Armazena os dados na struct
    data_struct.time = t_exp;
    data_struct.signals_clean = s_obs_by_line_clean;
    data_struct.signals = s_obs_by_line_noisy; % Usa a nova versão com ruído aditivo
    exp_data{exp_idx} = data_struct;
    %Plotting
    t_plot = exp_data{exp_idx}.time;
    s_plot_by_line = exp_data{exp_idx}.signals_clean; % Plot the clean data
    s_plot_global = sum(s_plot_by_line, 2);
    
    fig2 = figure('Position', [100, 100, 1200, 600]);
    subplot(1,2,1); plot(t_exp, sum(s_obs_by_line_clean, 2), 'r-', 'LineWidth', 2); title('Global (Summed) Sensorgram'); xlabel('Time (s)'); ylabel('Total s_{obs}(t)'); grid on;
    subplot(1,2,2); plot(t_exp, s_obs_by_line_clean, 'LineWidth', 1.5); title('Line-by-Line Sensorgrams (1D Heterogeneity)'); xlabel('Time (s)'); ylabel('s_{obs, j}(t)'); legend(arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false)); grid on;
    sgtitle(sprintf('Experiment %d: Ground Truth Data (RU)', exp_idx));
    
    % --- SAVE FIGURE ---
    print(fig2, sprintf('Adsorption/LineAverageModel/Figures/figure_2_ru_sensorgrams_exp%d.png', exp_idx), '-dpng', '-r300');
    print(fig2, sprintf('Adsorption/LineAverageModel/Figures/EPS/figure_2_ru_sensorgrams_exp%d.eps', exp_idx), '-depsc');
    % 1. Calculate total rows by summing up data points from all experiments
    total_rows = total_rows + numel(exp_data{exp_idx}.signals_clean);
end
fprintf('Data generation complete.\n');    
toc;
%% --- STEP 2: IDENTIFIABILITY ANALYSIS (1D HETEROGENEITY) ---
fprintf('\nStarting 1D Identifiability Analysis...\n');
% --- MODIFIED --- Create the "true" 1D parameter vector
% We assume the true line parameter is the average over the flow direction
kon_ads_2D = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
koff_ads_2D = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
smax_ads_2D = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
% Average parameters along the x-axis (flow direction) to get line parameters
p_true_kon_1D = mean(kon_ads_2D, 1)';   % [ads_y_dim x 1]
p_true_koff_1D = mean(koff_ads_2D, 1)'; % [ads_y_dim x 1]
p_true_smax_1D = mean(smax_ads_2D, 1)'; % [ads_y_dim x 1]
p_true_1D = [p_true_kon_1D; p_true_koff_1D; p_true_smax_1D];
N_params_1D = length(p_true_1D);
% --- START: PREALLOCATION FOR JACOBIAN ---
tic;
% 2. Preallocate the combined Jacobian matrix
J_combined = zeros(total_rows, N_params_1D);
% 3. Initialize a row indexer
current_row_start = 1;
% --- END: PREALLOCATION FOR JACOBIAN ---
% Loop to compute and fill the Jacobian
for exp_idx = 1:n_exp
    setting = exp_settings(exp_idx);
    fprintf('Computing Jacobian for Experiment %d (1D model)...\n', exp_idx);
    
    % REUSE baseline output from the data we already generated. No need to re-run.
    s_obs_matrix = exp_data{exp_idx}.signals_clean;
    s_obs_vector = s_obs_matrix(:); % Vectorize for Jacobian
    
    n_rows_exp = length(s_obs_vector);
    J_exp = zeros(n_rows_exp, N_params_1D);
    h = 1e-5;
    
    parfor k = 1:N_params_1D
        p_pert = p_true_1D;
        p_pert(k) = p_pert(k) * (1 + h);
        
        % We only need to run the perturbed simulation here
        [~, s_pert_matrix] = run_single_experiment_1D_model(p_pert, setting, model_config);
        s_pert_vector = s_pert_matrix(:);
        
        J_exp(:, k) = (s_pert_vector - s_obs_vector) / (p_true_1D(k) * h);
    end
    
    % --- Fill the preallocated matrix ---
    row_range = current_row_start : (current_row_start + n_rows_exp - 1);
    J_combined(row_range, :) = J_exp;
    current_row_start = current_row_start + n_rows_exp;
end
toc;
% ================= START OF SVD ANALYSIS =================
rankJ = rank(J_combined);
fprintf('\n1D Identifiability Analysis Results:\n');
fprintf('Total Parameters (3 * N_y): %d\n', N_params_1D);
fprintf('Rank of Combined Jacobian: %d\n', rankJ);
if rankJ < N_params_1D
    fprintf('WARNING: The model is structurally unidentifiable. Rank < Number of Parameters.\n');
else
    fprintf('SUCCESS: The model appears to be structurally identifiable (Jacobian has full rank).\n');
end
fprintf('Now performing SVD analysis to investigate practical identifiability...\n');
% --- Step 1: Perform Singular Value Decomposition ---
% Use the 'econ' flag for efficiency, as we only need the first N_params_1D vectors
[~, S, V] = svd(J_combined, 'econ');
% Extract the diagonal singular values
singular_values = diag(S);
% --- Step 2: Analyze and Plot Singular Values ---
svg_fig = figure('Name', 'SVD Analysis of Jacobian', 'Position', [100, 100, 1400, 600]);
subplot(1, 2, 1);
semilogy(singular_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Singular Values of the Jacobian');
xlabel('Singular Value Index');
ylabel('Magnitude (log scale)');
xlim([0, N_params_1D + 1]);
% Add text for the condition number
cond_number = singular_values(1) / singular_values(end);
legend(sprintf('Condition Number: %.2e', cond_number));
% --- Step 3: Analyze and Plot Parameter Combinations (Right Singular Vectors) ---
subplot(1, 2, 2);
imagesc(abs(V)); % Use absolute value for clarity of magnitude
colorbar;
title('Parameter Contributions to Singular Vectors (V)');
xlabel('Singular Vector Index (1=Most Identifiable -> N=Least Identifiable)');
ylabel('Parameter Index');
% Create meaningful labels for the y-axis
param_labels = [arrayfun(@(i) sprintf('kon_{%d}', i), 1:ads_y_dim, 'UniformOutput', false), ...
                arrayfun(@(i) sprintf('koff_{%d}', i), 1:ads_y_dim, 'UniformOutput', false), ...
                arrayfun(@(i) sprintf('smax_{%d}', i), 1:ads_y_dim, 'UniformOutput', false)];
yticks(1:N_params_1D);
yticklabels(param_labels);
%sgtitle('SVD-based Identifiability Analysis', 'FontSize', 16, 'FontWeight', 'bold');
print(svg_fig, 'Adsorption/LineAverageModel/Figures/figure_svd_analysis.png', '-dpng', '-r300');
print(svg_fig, 'Adsorption/LineAverageModel/Figures/EPS/figure_svd_analysis.eps', '-depsc');
% ================= END OF SVD ANALYSIS =================
%% --- STEP 3: PARAMETER IDENTIFICATION (1D HETEROGENEITY) ---   
fprintf('\nStarting 1D Parameter Identification...\n');
% --- MODIFIED --- Setup initial guess and bounds for 1D model
homog_kon = homog_params(1);
homog_koff = homog_params(2);
homog_smax_per_line = (homog_params(3) / ads_y_dim); % Evenly distributed total smax
p0_base_kon = log10(homog_kon) * ones(ads_y_dim, 1);
p0_base_koff = log10(homog_koff) * ones(ads_y_dim, 1);
p0_base_smax = log10(homog_smax_per_line) * ones(ads_y_dim, 1);
p0_base = [p0_base_kon; p0_base_koff; p0_base_smax];
    % Define TIGHT bounds in log10 space to keep the optimizer in a stable region.
% These correspond to kon=[1e2, 1e7], koff=[1e-5, 1], smax_per_line=[1e-4, 5]
lb_kon = log10(1e3);  ub_kon = log10(1e7);
lb_koff = log10(1e-6); ub_koff = log10(1e-1); 
lb_smax = log10(1e-4); ub_smax = log10(3000);   
lb = [repmat(lb_kon, ads_y_dim, 1); ...
      repmat(lb_koff, ads_y_dim, 1); ...
      repmat(lb_smax, ads_y_dim, 1)];
ub = [repmat(ub_kon, ads_y_dim, 1); ...
      repmat(ub_koff, ads_y_dim, 1); ...
      repmat(ub_smax, ads_y_dim, 1)];
% 3. Create the final p0 by adding a small random perturbation AND clamping to the bounds
rng('default'); % For reproducible randomness
noise_level = 0.1; % Small perturbation (e.g., 0.1 standard deviations in log space)
p0_rand = p0_base + noise_level * randn(size(p0_base));
% Clamp the randomized p0 to be within the bounds
p0 = max(lb, p0_rand); % Enforce lower bound
p0 = min(ub, p0);     % Enforce upper bound
% --- NEW --- Create a handle to the plotter function with all necessary data
optim_plot_fun = @(log_p, optim_v, state) optimPlotter_1D(...
    log_p, optim_v, state, ...
    p_true_1D, exp_data, ads_y_dim);
% --- MODIFIED --- Add the 'OutputFcn' to your optimization options
optim_opts = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'MaxIterations', 10, ...
    'UseParallel', true, ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'OutputFcn', optim_plot_fun); % This tells lsqnonlin to call our plotter
    
% --- MODIFIED --- The residual function remains the same
residual_fun = @(log_params) compute_residuals_1D_model(log_params, exp_settings, exp_data, model_config);
tic;
% Run optimization
[opt_log_params, resnorm, residual, exitflag, output, lambda, J_opt] = lsqnonlin(residual_fun, p0, lb, ub, optim_opts);
% --- NEW: STEP 3.5: CALCULATE CONFIDENCE INTERVALS ---
fprintf('\nCalculating 95%% Confidence Intervals for the estimated parameters...\n');

% 1. Estimate the variance of the measurement noise from the residuals
% Degrees of freedom = (num_data_points - num_parameters)
dof = numel(residual) - N_params_1D; 
noise_variance_est = resnorm / dof; % This is our estimate for sigma_s^2

% 2. Calculate the covariance matrix of the parameters (in log10 space)
% The formula is Cov(p) = sigma_s^2 * inv(J' * J)
% We use pinv (pseudo-inverse) for better numerical stability than inv
covariance_matrix_log = noise_variance_est * pinv(full(J_opt' * J_opt));

% 3. Extract the variances and standard errors for each parameter
param_variances_log = diag(covariance_matrix_log);
param_stderr_log = sqrt(param_variances_log);

% 4. Calculate the 95% confidence intervals
% For a 95% CI, the critical value is ~1.96
ci_95_log = 1.96 * param_stderr_log;

% The final results are the optimal log parameters +/- the interval
upper_bound_log = opt_log_params + ci_95_log;
lower_bound_log = opt_log_params - ci_95_log;

fprintf('Confidence interval calculation complete.\n');

toc;
% --- MODIFIED --- Analyze and plot results for 1D model
opt_params_1D = 10.^opt_log_params;
% Plot recovery of the 1D parameters
% Convert intervals back to linear space for plotting
upper_bound_lin = 10.^upper_bound_log;
lower_bound_lin = 10.^lower_bound_log;
y_errors_pos = upper_bound_lin - opt_params_1D;
y_errors_neg = opt_params_1D - lower_bound_lin;

% Update the call to the plotting function
plot_parameter_recovery_1D(p_true_1D, opt_params_1D, 10.^p0, ads_y_dim, y_errors_neg, y_errors_pos);% =========================================================================
%% --- STEP 4: PHYSICALLY-ACCURATE VALIDATION WORKFLOW ---
% =========================================================================
fprintf('\n--- Starting Full Physical Model Validation ---\n');
% --- Step 1: Define Optical and Physical Constants ---
fprintf('Defining optical parameters...\n');
wavelength = 670; % nm
d1 = 50;          % Gold film thickness (nm)
% The analyte layer (n2) thickness is effectively infinite for the evanescent wave
d2 = 1000;        % Effectively infinite analyte layer (nm)
n0 = sqrt(2.3104);         % Optical substrate (Prism)
n1 = sqrt(-14.379 + 1.0084j); % Gold film (complex RI)
n_bulk = sqrt(1.7876);         % Flow cell solution (baseline buffer)
% Define the angular range for SPR curve calculation
angle_range = linspace(65, 80, 1280); % [start_angle, end_angle, num_points]
% Define the conversion factor from Response Units (RU) to Refractive Index Units (RIU)
% 1000 RU = 0.001 RIU change
RU_TO_RIU = 0.001 / 1000;
% --- Step 2: Get the Ground-Truth Sensorgram Data (in RU) ---
% We use the clean, line-averaged data from our initial simulation
t_exp = exp_data{1}.time;
s_obs_ru = exp_data{1}.signals_clean; % Sensorgrams in RU
% --- Step 3: Convert Sensorgrams to Resonance Angles via Fresnel Model ---
fprintf('Processing %d time points for %d lines...\n', size(s_obs_ru, 1), size(s_obs_ru, 2));
%% 
% Preallocate matrices to store the calculated results
theta_spr_vs_time = zeros(size(s_obs_ru));
formula_response_vs_time = zeros(size(s_obs_ru));
baseline_offset = calculate_sensorgram_from_formula(n_bulk, n_bulk, n1, d2, wavelength);
fprintf('Calculated baseline offset of %.4f will be subtracted.\n', baseline_offset);
video_frames_folder = 'Adsorption/LineAverageModel/spr_video_frames';
% --- Setup for Video Frame Generation (only if flag is true) ---
if generate_video_frames
    fprintf('Generating SPR IMAGES FOR VIDEO...\n');
    if ~exist(video_frames_folder, 'dir')
        mkdir(video_frames_folder);
    else
        delete(fullfile(video_frames_folder, '*.png'));
    end
end
% --- Main Calculation and Optional Frame Generation Loop ---
tic;
num_time_points = size(s_obs_ru, 1);
fprintf('Processing %d time points...\n', num_time_points);
n2_vs_time = zeros(size(s_obs_ru)); 
for t_idx = 1:size(s_obs_ru, 1)
    spr_image_matrix = zeros(ads_y_dim, length(angle_range));
    % Loop through each line
   parfor j_idx = 1:size(s_obs_ru, 2)
        current_ru = s_obs_ru(t_idx, j_idx);
        n2_analyte = n_bulk + (current_ru * RU_TO_RIU);
        n2_vs_time(t_idx, j_idx) = n2_analyte; % Salva o RI calculado
        % Calculate full SPR curve, resonance angle, and formula response
        [Rp_curve, resonance_angle] = fresnel_spr_curve(angle_range, n0, n1, n2_analyte, n_bulk, d1, d2, wavelength);
        spr_image_matrix(j_idx, :) = Rp_curve;
        theta_spr_vs_time(t_idx, j_idx) = resonance_angle;
        delta_neff = calculate_sensorgram_from_formula(n2_analyte, n_bulk, n1, d2, wavelength);
        formula_response_vs_time(t_idx, j_idx) = delta_neff - baseline_offset;
    end
    % --- This block for creating and saving images ONLY runs if the flag is true ---
    if generate_video_frames
        % Update the data of the existing image object
        image_to_save_8bit = uint8(spr_image_matrix * 255);
        
        % Salva a matriz de dados diretamente como uma imagem PNG de 8-bit.
        filename = fullfile(video_frames_folder, sprintf('frame_%04d.png', t_idx));
        imwrite(image_to_save_8bit, filename);
    end
    
    % Display progress
    if mod(t_idx, 100) == 0
        fprintf('Processed frame %d de %d...\n', t_idx, num_time_points);
    end
end
toc;    
fprintf('Generating advanced SPR curve and image visualization for publication...\n');
    
% --- Define common plotting properties for publication ---
target_fig_width_cm = 8.4; % IEEE single column width
base_font_size = 8;        % Match your paper's caption font size (e.g., 8pt)
line_width = 1.5;          % Line width for main plot lines
line_width_thin = 1.0;     % Line width for annotations (like xline)

% --- Calculations (from your existing code) ---
line_to_plot = round(ads_y_dim/2);
[~, t_idx_max_response] = max(s_obs_ru(:, line_to_plot));
n2_baseline = n_bulk;
[Rp_baseline, theta_res_baseline] = fresnel_spr_curve(angle_range, n0, n1, n2_baseline, n_bulk, d1, d2, wavelength);
ru_max = s_obs_ru(t_idx_max_response, line_to_plot);
n2_analyte_max = n_bulk + (ru_max * RU_TO_RIU);
[Rp_analyte, theta_res_analyte] = fresnel_spr_curve(angle_range, n0, n1, n2_analyte_max, n_bulk, d1, d2, wavelength);
angle_shift = theta_res_analyte - theta_res_baseline;

% --- Create the Figure and Manually Position Axes ---
fig3_pub = figure('Name', 'SPR Curve Shift with Image Visualization');

pos_main_plot = [0.15, 0.35, 0.77, 0.55]; 
pos_img1 = [0.15, 0.22, 0.77, 0.06];      
pos_img2 = [0.15, 0.15, 0.77, 0.06];      

ax_main = axes('Position', pos_main_plot);
ax_img_baseline = axes('Position', pos_img1);
ax_img_analyte = axes('Position', pos_img2);

% --- Plot 1: The Main SPR Curves (on the top axes) ---
% CHANGE 1: Shorten the DisplayName text
plot(ax_main, angle_range, Rp_baseline, 'b-', 'LineWidth', line_width, 'DisplayName', 'Baseline');
hold(ax_main, 'on');
plot(ax_main, angle_range, Rp_analyte, 'r-', 'LineWidth', line_width, 'DisplayName', 'Max Response');
xline(ax_main, theta_res_baseline, 'b--', 'LineWidth', line_width_thin, 'HandleVisibility', 'off');
xline(ax_main, theta_res_analyte, 'r--', 'LineWidth', line_width_thin, 'HandleVisibility', 'off');
hold(ax_main, 'off');
grid(ax_main, 'on');
box(ax_main, 'on');
ylabel(ax_main, 'Reflectivity', 'FontSize', base_font_size);
ylim(ax_main, [0, 1]);
set(ax_main, 'XTickLabel', []); 
set(ax_main, 'FontSize', base_font_size - 1);

% CHANGE 2: Move the n2 information into the title
% title_str = sprintf('SPR Curve Shift (n_{2,base}=%.4f)', n2_baseline);
% title(ax_main, title_str, 'FontSize', base_font_size);

% CHANGE 3: Create a compact, boxless legend
lgd = legend(ax_main, 'Location', 'southeast', 'FontSize', base_font_size - 2);
lgd.Box = 'off'; % This removes the box around the legend

% --- Plot 2 & 3: The "SPR Images" (code is unchanged) ---
imagesc(ax_img_baseline, angle_range, 1, Rp_baseline);
colormap(ax_img_baseline, 'gray'); caxis(ax_img_baseline, [0,1]);
set(ax_img_baseline, 'YTick', []); set(ax_img_baseline, 'XTickLabel', []);
imagesc(ax_img_analyte, angle_range, 1, Rp_analyte);
colormap(ax_img_analyte, 'gray'); caxis(ax_img_analyte, [0,1]);
set(ax_img_analyte, 'YTick', []);
xlabel(ax_img_analyte, 'Incident Angle (degrees)', 'FontSize', base_font_size);
set(ax_img_analyte, 'FontSize', base_font_size - 1);

% --- Link all X-Axes together (code is unchanged) ---
linkaxes([ax_main, ax_img_baseline, ax_img_analyte], 'x');
xlim(ax_main, [angle_range(1), angle_range(end)]); 

% --- Annotations (code is unchanged) ---
ax_pos = get(ax_main, 'Position');
xlims = get(ax_main, 'XLim'); ylims = get(ax_main, 'YLim');
y_arrow = 0.5; p1_data = [theta_res_baseline, y_arrow]; p2_data = [theta_res_analyte, y_arrow];
x_arrow_norm = ( [p1_data(1), p2_data(1)] - xlims(1) ) / diff(xlims);
y_arrow_norm = ( [p1_data(2), p2_data(2)] - ylims(1) ) / diff(ylims);
x_arrow_fig = ax_pos(1) + x_arrow_norm * ax_pos(3);
y_arrow_fig = ax_pos(2) + y_arrow_norm * ax_pos(4);
annotation('doublearrow', x_arrow_fig, y_arrow_fig, 'LineWidth', line_width_thin, 'Color', 'k', 'HeadStyle', 'vback2', 'HeadSize', 6);
text_str = sprintf('\\Delta\\theta_{SPR} = %.3f°', angle_shift);
text(ax_main, xlims(1) + 0.05*diff(xlims), ylims(1) + 0.9*diff(ylims), text_str, 'FontSize', base_font_size - 1, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'VerticalAlignment', 'top');

% --- SAVE THE FINAL FIGURE ---
save_pub_fig(fig3_pub, 'Adsorption/LineAverageModel/Figures/figure_spr_shift_composite_pub', target_fig_width_cm);
close(fig3_pub);
%% --- STEP 5: VIDEO CREATION ---
% =========================================================================
if generate_video_frames
    fprintf('\n--- Iniciando a criação do vídeo a partir dos frames salvos ---\n');
    tic;
    % --- THE FIX IS HERE (Part 1): Change the video profile and filename ---
    video_filename = 'Adsorption/LineAverageModel/Videos/spri_simulation_final.avi';
    outputVideo = VideoWriter(video_filename, 'Motion JPEG AVI');
    
    outputVideo.FrameRate = 0.5;
    outputVideo.Quality = 100; % Quality for AVI is 0-100
    open(outputVideo);
    % --- Pega e ordena a lista de todos os arquivos de imagem ---
    image_files_struct = dir(fullfile(video_frames_folder, '*.png'));
    image_files_cell = {image_files_struct.name};
    str_nums = regexp(image_files_cell, '\d+', 'match', 'once');
    num_vals = str2double(str_nums);
    [~, sorted_indices] = sort(num_vals);
    sorted_image_files = image_files_cell(sorted_indices);
    
    fprintf('Lendo %d frames para criar o vídeo...\n', length(sorted_image_files));
    for i = 1:length(sorted_image_files)
        img_path = fullfile(video_frames_folder, sorted_image_files{i});
        img = imread(img_path);
        
        % --- THE FIX IS HERE (Part 2): Ensure image is 3-channel RGB ---
        % The Motion JPEG AVI codec expects a 3-channel (RGB) frame.
        % If the image is grayscale, we replicate it into 3 channels.
        if size(img, 3) == 1
            img = cat(3, img, img, img);
        end
        
        % Escreve o frame no vídeo
        writeVideo(outputVideo, img);
    end
    % --- Finaliza e fecha o arquivo de vídeo ---
    close(outputVideo);
    toc;
    fprintf('\nVídeo salvo com sucesso como "%s".\n', video_filename);
end
%% --- STEP 6: PLOTS ---
% We plot the Resonance Angle directly. To validate, we overlay
% the original RU data on a second y-axis to show the shapes match.
%% 
fig_spr_curves = figure('Name', 'SPR Curves for All Lines at Max Response');
hold on;
for j_idx = 1:ads_y_dim
    plot(angle_range, spr_image_matrix(j_idx, :), 'LineWidth', 2, 'DisplayName', sprintf('Linha %d', j_idx));
end
hold off;
grid on;
xlabel('Ângulo de Incidência (graus)');
ylabel('Refletividade');
legend('Location', 'best');
ylim([0, 1]);
% --- SAVE FIGURE ---
print(fig_spr_curves, 'Adsorption/LineAverageModel/Figures/figure_spr_curves_all_lines.png', '-dpng', '-r300');
print(fig_spr_curves, 'Adsorption/LineAverageModel/Figures/EPS/figure_spr_curves_all_lines.eps', '-depsc');
%% 
fprintf('Generating publication-ready sensorgram comparison plot...\n');

% --- Define common plotting properties for publication ---
target_fig_width_cm = 8.4; % IEEE single column width
base_font_size = 8;        % Match your paper's caption font size (e.g., 8pt)
line_width = 1.2;          % Line width for plot lines

% --- Create the figure and layout ---
fig3b_pub = figure('Name', 'Proper Sensorgram: Resonance Angle vs. Time');
lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);

% Use tiledlayout for compact spacing
t = tiledlayout(1, length(lines_to_plot), 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Loop through the lines to plot ---
for i = 1:length(lines_to_plot)
    ax = nexttile; % Get the axes for the current tile
    line_idx = lines_to_plot(i);
    
    % --- Plot on the LEFT Y-axis (Resonance Angle) ---
    yyaxis(ax, 'left');
    plot(ax, t_exp, theta_spr_vs_time(:, line_idx), 'r-', 'LineWidth', line_width, 'DisplayName', 'Resonance Angle');
    
    hold(ax, 'on');
    grid(ax, 'on');
    
    % Style the left axis
    ax.YAxis(1).Color = 'r';
    if i == 1 % Only show the first y-label
        ylabel(ax, 'Res. Angle (deg)', 'FontSize', base_font_size);
    else
        % Remove redundant tick labels from interior plots
        set(ax, 'YTickLabel', []); 
    end
    
    % --- Plot on the RIGHT Y-axis (Response Units) ---
    yyaxis(ax, 'right');
    plot(ax, t_exp, s_obs_ru(:, line_idx), 'b--', 'LineWidth', line_width, 'DisplayName', 'Original (RU)');
    
    hold(ax, 'off');
    
    % Style the right axis
    ax.YAxis(2).Color = 'b';
    if i == length(lines_to_plot) % Only show the last y-label
        ylabel(ax, 'Response (RU)', 'FontSize', base_font_size);
    else
        set(ax, 'YTickLabel', []);
    end
    
    % --- Common properties for this subplot ---
    title(ax, sprintf('Line %d', line_idx), 'FontSize', base_font_size);
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1); % Set tick font size
    box(ax, 'on');
end

% --- Add shared elements to the entire layout ---
xlabel(t, 'Time (s)', 'FontSize', base_font_size);
lgd = legend(ax, 'Location', 'best', 'FontSize', base_font_size - 2);
lgd.Box = 'off';

%title(t, 'Final Validation: Physically Correct Sensorgram', 'FontSize', base_font_size + 1);
% --- SAVE THE FINAL FIGURE ---
save_pub_fig(fig3b_pub, 'Adsorption/LineAverageModel/Figures/figure_angle_vs_ru_pub', target_fig_width_cm);
close(fig3b_pub);
%% 
fprintf('Generating publication-ready final sensorgram plot...\n');

% --- Define common plotting properties for publication ---
target_fig_width_cm = 8.4; % IEEE single column width
base_font_size = 8;        % Match your paper's caption font size (e.g., 8pt)
line_width_thin = 0.8;     % Thinner lines for the dense plot
line_width_thick = 1.5;    % Thicker line for the single plot

% --- Calculate data ---
global_sensorgram = sum(formula_response_vs_time, 2);

% --- Create the Figure and Layout ---
fig_formula_pub = figure('Name', 'Final Sensorgram Results from Formula');
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Plot 1: All Individual Line Sensorgrams ---
ax1 = nexttile;
plot(ax1, t_exp, formula_response_vs_time, 'LineWidth', line_width_thin);
grid(ax1, 'on');
box(ax1, 'on');
title(ax1, 'Individual Lines', 'FontSize', base_font_size);
xlabel(ax1, 'Time (s)', 'FontSize', base_font_size);
ylabel(ax1, 'Change from Baseline (\DeltaN_s^{eff})', 'FontSize', base_font_size);
xlim(ax1, [0, t_exp(end)]);
set(ax1, 'FontSize', base_font_size - 1);

% Optional: Add a legend if you have a small number of lines
if ads_y_dim <= 5 % Reduced threshold for a compact plot
    legend(ax1, arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), ...
           'Location', 'northwest', 'FontSize', base_font_size - 2);
end

% --- Plot 2: Global (Summed) Sensorgram ---
ax2 = nexttile;
plot(ax2, global_sensorgram, 'r-', 'LineWidth', line_width_thick);
grid(ax2, 'on');
box(ax2, 'on');
title(ax2, 'Global (Summed)', 'FontSize', base_font_size);
xlabel(ax2, 'Time (s)', 'FontSize', base_font_size);
ylabel(ax2, 'Total Change (\Sigma\DeltaN_s^{eff})', 'FontSize', base_font_size);
xlim(ax2, [0, t_exp(end)]);
set(ax2, 'FontSize', base_font_size - 1);

% --- Add a main title to the layout ---
title(t, 'Sensorgrams from Analytical Formula', 'FontSize', base_font_size + 1);

% --- SAVE THE FINAL FIGURE ---
save_pub_fig(fig_formula_pub, 'Adsorption/LineAverageModel/Figures/figure_final_sensorgrams_formula_pub', target_fig_width_cm);
close(fig_formula_pub);

fprintf('Generating publication-ready absolute RI sensorgram plot...\n');

% --- Define common plotting properties for publication ---
target_fig_width_cm = 8.4; % IEEE single column width
base_font_size = 8;        % Match your paper's caption font size (e.g., 8pt)
line_width_thin = 0.8;     % Thinner lines for the dense plot
line_width_thick = 1.5;    % Thicker line for the single plot

% --- Step 1: Calculate data ---
absolute_neff_vs_time = formula_response_vs_time + n_bulk;
global_absolute_neff = n_bulk + global_sensorgram;

% --- Step 2: Create the Figure and Layout ---
fig_abs_ri_pub = figure('Name', 'Absolute Effective RI Sensorgrams');
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Plot 1: All Individual Line Sensorgrams ---
ax1 = nexttile;
plot(ax1, t_exp, absolute_neff_vs_time, 'LineWidth', line_width_thin);
grid(ax1, 'on');
box(ax1, 'on');
title(ax1, 'Individual Lines', 'FontSize', base_font_size);
xlabel(ax1, 'Time (s)', 'FontSize', base_font_size);
ylabel(ax1, 'Abs. Eff. RI (N_s^{eff})', 'FontSize', base_font_size);
xlim(ax1, [0, t_exp(end)]);
set(ax1, 'FontSize', base_font_size - 1);

% Optional: Add a legend if you have a small number of lines
if ads_y_dim <= 5 % Reduced threshold for a compact plot
    legend(ax1, arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), ...
           'Location', 'northwest', 'FontSize', base_font_size - 2);
end

% --- Plot 2: Global (Summed) Sensorgram ---
ax2 = nexttile;
plot(ax2, global_absolute_neff, 'r-', 'LineWidth', line_width_thick);
grid(ax2, 'on');
box(ax2, 'on');
title(ax2, 'Global (Summed)', 'FontSize', base_font_size);
xlabel(ax2, 'Time (s)', 'FontSize', base_font_size);
ylabel(ax2, 'Abs. Eff. RI (\Sigma N_s^{eff})', 'FontSize', base_font_size);
xlim(ax2, [0, t_exp(end)]);
set(ax2, 'FontSize', base_font_size - 1);

% --- Add a main title to the layout ---
title(t, 'Sensorgrams as Absolute Effective RI', 'FontSize', base_font_size + 1);

% --- SAVE THE FINAL FIGURE ---
save_pub_fig(fig_abs_ri_pub, 'Adsorption/LineAverageModel/Figures/figure_absolute_ri_sensorgrams_pub', target_fig_width_cm);
close(fig_abs_ri_pub);
% =========================================================================
%% --- STEP 7: INVERSE PROCESS (Corrected with Single Legend)---
% =========================================================================
% --- Step 7.1-7.3: Data Calculation (code is unchanged) ---
theta_extracted_from_frames = analyze_spr_frames_to_get_sensorgram(video_frames_folder, angle_range);
fprintf('Converting extracted angles back to Refractive Index via interpolation...\n');
n2_reconstructed_vs_time = zeros(size(theta_extracted_from_frames));
parfor j_idx = 1:ads_y_dim
    [unique_thetas, unique_indices] = unique(theta_spr_vs_time(:, j_idx));
    unique_n2s = n2_vs_time(unique_indices, j_idx);
    n2_reconstructed_vs_time(:, j_idx) = interp1(unique_thetas, unique_n2s, theta_extracted_from_frames(:, j_idx), 'linear', 'extrap');
end
ru_reconstructed = (n2_reconstructed_vs_time - n_bulk) / RU_TO_RIU;

% =====================================================================
% --- FINAL VALIDATION PLOTS WITH SINGLE, CLEAN LEGEND ---
% =====================================================================
lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);
target_fig_width_cm = 8.4; 
base_font_size = 8;        
line_width_thick = 1.5;    
line_width_thin = 1.0;     

%% --- Plot A: Validation in Resonance Angle units (Corrected) ---
fprintf('\n--- Angle Validation Results ---\n');
fig4a_pub = figure('Name', 'Final Validation: Original vs. Extracted Angle');
t_a = tiledlayout(1, length(lines_to_plot), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:length(lines_to_plot)
    ax = nexttile;
    line_idx = lines_to_plot(i);
    original_angle = theta_spr_vs_time(:, line_idx);
    extracted_angle = theta_extracted_from_frames(:, line_idx);
    
    % --- FIX: Print errors to the command window ---
    [offset, error_abs] = calculate_validation_error(original_angle, extracted_angle);
    fprintf('Angle Validation for Line %d: Abs. Error = %.2e (deg)\n', line_idx, error_abs);
    
    % --- FIX: Use static legend text ---
    plot(ax, t_exp, original_angle - original_angle(1), 'b-', 'LineWidth', line_width_thick, 'DisplayName', 'Original');
    hold(ax, 'on');
    plot(ax, t_exp, extracted_angle - extracted_angle(1), 'r--', 'LineWidth', line_width_thin, 'DisplayName', 'Reconstructed');
    hold(ax, 'off');
    
    grid(ax, 'on'); box(ax, 'on');
    title(ax, sprintf('Line %d', line_idx), 'FontSize', base_font_size);
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1);
    if i == 1, ylabel(ax, '\Delta Angle (deg)', 'FontSize', base_font_size); end
end
% --- FIX: Create a single, shared legend for the entire figure ---
lgd = legend(ax); % Create legend attached to the LAST axes
lgd.Layout.Tile = 'east'; % Move the legend to its own space outside the plots

% title(t_a, 'Validation in Resonance Angle Space', 'FontSize', base_font_size + 1);
xlabel(t_a, 'Time (s)', 'FontSize', base_font_size);
save_pub_fig(fig4a_pub, 'Adsorption/LineAverageModel/Figures/figure_4a_validation_angle_pub', target_fig_width_cm);
close(fig4a_pub);

%% --- Plot B: Validation in Refractive Index units (Corrected) ---
fprintf('\n--- Refractive Index Validation Results ---\n');
fig4b_pub = figure('Name', 'Final Validation: Original vs. Reconstructed RI');
t_b = tiledlayout(1, length(lines_to_plot), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:length(lines_to_plot)
    ax = nexttile;
    line_idx = lines_to_plot(i);
    original_ri = n2_vs_time(:, line_idx);
    reconstructed_ri = n2_reconstructed_vs_time(:, line_idx);
    
    % --- FIX: Print errors to the command window ---
    [offset, error_abs] = calculate_validation_error(original_ri, reconstructed_ri);
    fprintf('RI Validation for Line %d: Abs. Error = %.2e (RIU)\n', line_idx, error_abs);

    % --- FIX: Use static legend text ---
    plot(ax, t_exp, original_ri - original_ri(1), 'b-', 'LineWidth', line_width_thick, 'DisplayName', 'Original');
    hold(ax, 'on');
    plot(ax, t_exp, reconstructed_ri - reconstructed_ri(1), 'r--', 'LineWidth', line_width_thin, 'DisplayName', 'Reconstructed');
    hold(ax, 'off');
    
    grid(ax, 'on'); box(ax, 'on');
    title(ax, sprintf('Line %d', line_idx), 'FontSize', base_font_size);
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1);
    if i == 1, ylabel(ax, '\DeltaRIU', 'FontSize', base_font_size); end
end
% --- FIX: Create a single, shared legend for the entire figure ---
lgd = legend(ax);
lgd.Layout.Tile = 'east';

% title(t_b, 'Validation in Refractive Index Space', 'FontSize', base_font_size + 1);
xlabel(t_b, 'Time (s)', 'FontSize', base_font_size);
save_pub_fig(fig4b_pub, 'Adsorption/LineAverageModel/Figures/figure_4b_validation_ri_pub', target_fig_width_cm);
close(fig4b_pub);

%% --- Plot C: Validation in Response Units (Corrected) ---
fprintf('\n--- Response Unit Validation Results ---\n');
fig4c_pub = figure('Name', 'Final Round-Trip Validation: RU Original vs. Reconstructed');
t_c = tiledlayout(1, length(lines_to_plot), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:length(lines_to_plot)
    ax = nexttile;
    line_idx = lines_to_plot(i);
    original_ru = s_obs_ru(:, line_idx);
    reconstructed_ru_line = ru_reconstructed(:, line_idx);

    % --- FIX: Print errors to the command window ---
    [offset, error_abs] = calculate_validation_error(original_ru, reconstructed_ru_line);
    fprintf('RU Validation for Line %d: Abs. Error = %.2e (RU)\n', line_idx, error_abs);

    % --- FIX: Use static legend text ---
    plot(ax, t_exp, original_ru - original_ru(1), 'b-', 'LineWidth', line_width_thick, 'DisplayName', 'Original');
    hold(ax, 'on');
    plot(ax, t_exp, reconstructed_ru_line - reconstructed_ru_line(1), 'r--', 'LineWidth', line_width_thin, 'DisplayName', 'Reconstructed');
    hold(ax, 'off');

    grid(ax, 'on'); box(ax, 'on');
    title(ax, sprintf('Line %d', line_idx), 'FontSize', base_font_size);
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1);
    if i == 1, ylabel(ax, 'Response Units (RU)', 'FontSize', base_font_size); end
end
% --- FIX: Create a single, shared legend for the entire figure ---
lgd = legend(ax);
lgd.Layout.Tile = 'east';

% title(t_c, 'Final Round-Trip Validation in Response Units', 'FontSize', base_font_size + 1);
xlabel(t_c, 'Time (s)', 'FontSize', base_font_size);
save_pub_fig(fig4c_pub, 'Adsorption/LineAverageModel/Figures/figure_4c_validation_ru_pub', target_fig_width_cm);
close(fig4c_pub);


%% --- FINAL COMPOSITE FIGURE: Uniting All Validation Plots (3x3 Grid) ---
% =========================================================================
% This section creates a new, single figure by stacking the three 1x3 
% validation plots (Angle, RI, RU) vertically.
% =========================================================================
fprintf('\n\n--- Generating Final 3x3 Composite Validation Figure ---\n');

% --- Define Publication Style Parameters ---
% A 3x3 figure needs a wider format to be readable. 18cm is a common full-page width.
target_fig_width_cm = 18;  
base_font_size = 8;        
line_width_thick = 1.2;    
line_width_thin = 1.0;     
lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);

% --- Create the master figure and 3x3 layout ---
fig_summary = figure('Name', 'Comprehensive Validation Summary');
t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% =====================================================================
% --- ROW 1: Angle Validation ---
% =====================================================================
for i = 1:length(lines_to_plot)
    ax = nexttile;
    line_idx = lines_to_plot(i);
    original_angle = theta_spr_vs_time(:, line_idx);
    extracted_angle = theta_extracted_from_frames(:, line_idx);
    
    plot(ax, t_exp, original_angle - original_angle(1), 'b-', 'LineWidth', line_width_thick, 'DisplayName', 'Original');
    hold(ax, 'on');
    plot(ax, t_exp, extracted_angle - extracted_angle(1), 'r--', 'LineWidth', line_width_thin, 'DisplayName', 'Reconstructed');
    hold(ax, 'off');
    
    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1, 'XTickLabel', []); % Remove x-ticks from top row
    
    title(ax, sprintf('Line %d', line_idx), 'FontSize', base_font_size);
    if i == 1, ylabel(ax, '\Delta Angle (deg)', 'FontSize', base_font_size); end
end

% =====================================================================
% --- ROW 2: RI Validation ---
% =====================================================================
for i = 1:length(lines_to_plot)
    ax = nexttile;
    line_idx = lines_to_plot(i);
    original_ri = n2_vs_time(:, line_idx);
    reconstructed_ri = n2_reconstructed_vs_time(:, line_idx);

    plot(ax, t_exp, original_ri - original_ri(1), 'b-', 'LineWidth', line_width_thick, 'DisplayName', 'Original');
    hold(ax, 'on');
    plot(ax, t_exp, reconstructed_ri - reconstructed_ri(1), 'r--', 'LineWidth', line_width_thin, 'DisplayName', 'Reconstructed');
    hold(ax, 'off');
    
    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1, 'XTickLabel', []); % Remove x-ticks from middle row
    
    if i == 1, ylabel(ax, '\DeltaRIU', 'FontSize', base_font_size); end
end

% =====================================================================
% --- ROW 3: RU Validation ---
% =====================================================================
for i = 1:length(lines_to_plot)
    ax = nexttile;
    line_idx = lines_to_plot(i);
    original_ru = s_obs_ru(:, line_idx);
    reconstructed_ru_line = ru_reconstructed(:, line_idx);

    plot(ax, t_exp, original_ru - original_ru(1), 'b-', 'LineWidth', line_width_thick, 'DisplayName', 'Original');
    hold(ax, 'on');
    plot(ax, t_exp, reconstructed_ru_line - reconstructed_ru_line(1), 'r--', 'LineWidth', line_width_thin, 'DisplayName', 'Reconstructed');
    hold(ax, 'off');

    grid(ax, 'on'); box(ax, 'on');
    xlim(ax, [0, t_exp(end)]);
    set(ax, 'FontSize', base_font_size - 1); % Keep x-ticks on bottom row
    
    if i == 1, ylabel(ax, 'Response Units (RU)', 'FontSize', base_font_size); end
end

% --- Add Shared Legend and Labels ---
lgd = legend(ax); % Create legend attached to the LAST axes handle
lgd.Layout.Tile = 'South'; % Move the legend below all plots
lgd.NumColumns = 2;        % Arrange legend items horizontally
lgd.FontSize = base_font_size;

% Add a single, shared X-axis label to the whole layout
xlabel(t, 'Time (s)', 'FontSize', base_font_size + 1);

% Add a main title for the entire figure
% title(t, 'Comprehensive Round-Trip Validation', 'FontSize', base_font_size + 2, 'FontWeight', 'bold');

% --- Save the Final Composite Figure ---
save_pub_fig(fig_summary, 'Adsorption/LineAverageModel/Figures/figure_VALIDATION_summary_composite', target_fig_width_cm);
close(fig_summary);
% =========================================================================
% --- FINAL PUBLICATION PLOT: ANNOTATED SENSORGRAM FIT ---
% =========================================================================
fprintf('\nGenerating annotated experiment plot for publication...\n');

% --- Define common plotting properties for publication ---
target_fig_width_cm = 8.4; % IEEE single column width
base_font_size = 8;        % Match your paper's caption font size (e.g., 8pt)
line_width_fit = 1.5;      % Line width for the main fit
line_width_anno = 0.8;     % Line width for annotation lines
marker_size = 4;           % Size of the data markers

% --- Select and calculate data ---
exp_to_plot = 1;
line_to_plot = round(ads_y_dim/2);
setting_to_plot = exp_settings(exp_to_plot);
t_exp = exp_data{exp_to_plot}.time;
s_obs_noisy = exp_data{exp_to_plot}.signals(:, line_to_plot);
[~, s_final_fit_all_lines] = run_single_experiment_1D_model(opt_params_1D, setting_to_plot, model_config);
s_final_fit = s_final_fit_all_lines(:, line_to_plot);

% --- Create the Plot ---
fig_annotated_pub = figure('Name', 'Annotated Sensorgram Fit');
ax = gca; % Get current axes
hold(ax, 'on');

% --- Plot the noisy data and the final smooth fit ---
plot(ax, t_exp, s_obs_noisy, '.', 'Color', [0.6 0.6 1], 'MarkerSize', marker_size, 'DisplayName', 'Noisy Data');
plot(ax, t_exp, s_final_fit, 'r-', 'LineWidth', line_width_fit, 'DisplayName', 'Model Fit');

% --- Add vertical lines and text annotations for each phase ---
y_lims = ylim; % Get current y-axis limits
text_y_pos = y_lims(1) + 0.95 * (y_lims(2) - y_lims(1)); % Position text at 95% of y-axis
all_t_breaks = [0, setting_to_plot.pulse_times, t_total];
all_concs = [setting_to_plot.pulse_concs, c_diss];

for i = 1:length(all_concs)
    t_start = all_t_breaks(i);
    t_end = all_t_breaks(i+1);
    
    % Draw vertical line at the start of the new phase
    if i > 1
        line(ax, [t_start, t_start], y_lims, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', line_width_anno, 'HandleVisibility', 'off');
    end
    
    % Add text annotation in the middle of the phase
    text_x_pos = (t_start + t_end) / 2;
    text_str = sprintf('C=%.1eM', all_concs(i)); % Shortened text
    
    text(ax, text_x_pos, text_y_pos, text_str, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', base_font_size - 1, ... % Use a slightly smaller font for annotations
        'FontWeight', 'normal', ...          % Normal weight is better for small fonts
        'BackgroundColor', [1 1 1 0.7], ... 
        'Margin', 2); % Smaller margin for smaller text
end

% --- Finalize Plot ---
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
axis(ax, 'tight');
%title(ax, sprintf('Model Fit for Line %d', line_to_plot), 'FontSize', base_font_size);
xlabel(ax, 'Time (s)', 'FontSize', base_font_size);
ylabel(ax, 'Response (RU)', 'FontSize', base_font_size);
legend(ax, 'Location', 'southeast', 'FontSize', base_font_size - 1);
set(ax, 'FontSize', base_font_size - 1);

% --- SAVE THE FINAL FIGURE ---
save_pub_fig(fig_annotated_pub, 'Adsorption/LineAverageModel/Figures/figure_6_annotated_fit_pub', target_fig_width_cm);
close(fig_annotated_pub);



