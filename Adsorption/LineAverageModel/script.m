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
% Extract the 2D slice of each parameter for the entire grid
kon_slice = squeeze(kon_grid_heterog(:, :, ads_layer));
koff_slice = squeeze(koff_grid_heterog(:, :, ads_layer));
smax_slice = squeeze(smax_grid_heterog(:, :, ads_layer));
% Create the figure and subplots
fig1 = figure('Name', 'Surface Parameter Heatmaps', 'Position', [100, 100, 1600, 450]);
% --- Define the rectangle's position and size from your range variables ---
% Position format: [x_start, y_start, width, height]
% We subtract 0.5 to center the rectangle around the pixels.
rect_pos = [ads_x_range(1)-0.5, ads_y_range(1)-0.5, ...
            ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1];
            
% --- Define the explicit axes for the plot ---
x_axis_coords = 1:gridN_x;
y_axis_coords = 1:gridN_y;

% --- 2D Plot for k_on ---
ax1 = subplot(1, 3, 1);
imagesc(x_axis_coords, y_axis_coords, kon_slice'); % Use imagesc and transpose (') for intuitive orientation
axis xy; % Place the y-axis origin at the bottom-left
hold on; % Prepare to draw on top of the image
rectangle('Position', rect_pos, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
hold off;
title('k_{on} Surface');
xlabel('Position along Flow (x)');
ylabel('Line Index (y)');
colorbar;

% --- 2D Plot for k_off ---
ax2 = subplot(1, 3, 2);
imagesc(x_axis_coords, y_axis_coords, koff_slice');
axis xy;
hold on;
rectangle('Position', rect_pos, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
hold off;
title('k_{off} Surface');
xlabel('Position along Flow (x)');
ylabel('Line Index (y)');
colorbar;

% --- 2D Plot for s_max ---
ax3 = subplot(1, 3, 3);
imagesc(x_axis_coords, y_axis_coords, smax_slice');
axis xy;
hold on;
rectangle('Position', rect_pos, 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
hold off;
title('s_{max} Surface');
xlabel('Position along Flow (x)');
ylabel('Line Index (y)');
colorbar;

%sgtitle('2D Visualization of Surface Parameter Heterogeneity', 'FontSize', 16, 'FontWeight', 'bold');

% --- SAVE FULL FIGURE ---
print(fig1, 'Adsorption/LineAverageModel/Figures/figure_1_surf_params_2D.png', '-dpng', '-r300');
print(fig1, 'Adsorption/LineAverageModel/Figures/EPS/figure_1_surf_params_2D.eps', '-depsc');

% --- SAVE EACH SUBPLOT INDIVIDUALLY ---
fprintf('Saving individual surface parameter heatmaps as EPS files...\n');
base_path = 'Adsorption/LineAverageModel/Figures/EPS/';
save_subplot_as_eps(ax1, [base_path, 'surf_params_heatmap_kon.eps']);
save_subplot_as_eps(ax2, [base_path, 'surf_params_heatmap_koff.eps']);
save_subplot_as_eps(ax3, [base_path, 'surf_params_heatmap_smax.eps']);
fprintf('Finished saving individual surface parameter heatmaps.\n');
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
[opt_log_params, ~] = lsqnonlin(residual_fun, p0, lb, ub, optim_opts);
toc;
% --- MODIFIED --- Analyze and plot results for 1D model
opt_params_1D = 10.^opt_log_params;
% Plot recovery of the 1D parameters
plot_parameter_recovery_1D(p_true_1D, opt_params_1D, 10.^p0, ads_y_dim);
% =========================================================================
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
angle_range = linspace(65, 80, 1000); % [start_angle, end_angle, num_points]
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
fprintf('Generating advanced SPR curve and image visualization...\n');
    
% --- Calculations (from your existing code) ---
line_to_plot = round(ads_y_dim/2);
[~, t_idx_max_response] = max(s_obs_ru(:, line_to_plot));
n2_baseline = n_bulk;
[Rp_baseline, theta_res_baseline] = fresnel_spr_curve(angle_range, n0, n1, n2_baseline, n_bulk, d1, d2, wavelength);
ru_max = s_obs_ru(t_idx_max_response, line_to_plot);
n2_analyte_max = n_bulk + (ru_max * RU_TO_RIU);
[Rp_analyte, theta_res_analyte] = fresnel_spr_curve(angle_range, n0, n1, n2_analyte_max, n_bulk, d1, d2, wavelength);
% Calculate the shift values for annotation
angle_shift = theta_res_analyte - theta_res_baseline;
ri_change = n2_analyte_max - n2_baseline;
% --- Create the Figure and Manually Position Axes ---
fig3 = figure('Name', 'SPR Curve Shift with Image Visualization', 'Position', [100, 100, 800, 800]);
% Define positions for the plots: [left, bottom, width, height]
pos_main_plot = [0.13, 0.35, 0.77, 0.55]; % Large plot on top
pos_img1 = [0.13, 0.20, 0.77, 0.05];      % Thin image strip below
pos_img2 = [0.13, 0.10, 0.77, 0.05];      % Second thin image strip
% Create the axes objects
ax_main = axes('Position', pos_main_plot);
ax_img_baseline = axes('Position', pos_img1);
ax_img_analyte = axes('Position', pos_img2);
% --- Plot 1: The Main SPR Curves (on the top axes) ---
plot(ax_main, angle_range, Rp_baseline, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Baseline (n_2 = %.4f)', n2_baseline));
hold(ax_main, 'on');
plot(ax_main, angle_range, Rp_analyte, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Max Response (n_2 = %.4f)', n2_analyte_max));
xline(ax_main, theta_res_baseline, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(ax_main, theta_res_analyte, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold(ax_main, 'off');
grid(ax_main, 'on');
ylabel(ax_main, 'Reflectivity');
title(ax_main, sprintf('SPR Curve Shift for Line %d', line_to_plot));
legend(ax_main, 'Location', 'northeast');
ylim(ax_main, [0, 1]);
set(ax_main, 'XTickLabel', []); % Remove x-axis labels to avoid overlap
% --- Plot 2: The "SPR Image" for the Baseline Curve ---
imagesc(ax_img_baseline, angle_range, 1, Rp_baseline);
colormap(ax_img_baseline, 'gray');
caxis(ax_img_baseline, [0,1]);
set(ax_img_baseline, 'YTick', []); % Remove y-axis ticks
set(ax_img_baseline, 'XTickLabel', []); % Remove x-axis labels
% --- Plot 3: The "SPR Image" for the Max Response Curve ---
imagesc(ax_img_analyte, angle_range, 1, Rp_analyte);
colormap(ax_img_analyte, 'gray');
caxis(ax_img_analyte, [0,1]);
set(ax_img_analyte, 'YTick', []);
xlabel(ax_img_analyte, 'Incident Angle (degrees)'); % Only show x-label on the bottom plot
% --- Link all X-Axes together ---
linkaxes([ax_main, ax_img_baseline, ax_img_analyte], 'x');
xlim(ax_main, [angle_range(1), angle_range(end)]); % Set initial limits
% --- Add Annotations for the Shift Arrow and Text ---
% Get position of the main plot to convert data coordinates to figure coordinates
ax_pos = get(ax_main, 'Position');
xlims = get(ax_main, 'XLim');
ylims = get(ax_main, 'YLim');
% Arrow coordinates in data space
y_arrow = 0.5; % Y position for the arrow
p1_data = [theta_res_baseline, y_arrow];
p2_data = [theta_res_analyte, y_arrow];
% Convert to normalized figure units for the annotation
x_arrow_norm = ( [p1_data(1), p2_data(1)] - xlims(1) ) / diff(xlims);
y_arrow_norm = ( [p1_data(2), p2_data(2)] - ylims(1) ) / diff(ylims);
x_arrow_fig = ax_pos(1) + x_arrow_norm * ax_pos(3);
y_arrow_fig = ax_pos(2) + y_arrow_norm * ax_pos(4);
% Draw the double-headed arrow
annotation('doublearrow', x_arrow_fig, y_arrow_fig, 'LineWidth', 2, 'Color', 'k', 'HeadStyle', 'vback2', 'HeadSize', 10);
% Add text annotation with shift information
text_str = sprintf('\\Delta\\theta_{SPR} = %.3f°', angle_shift);
text(ax_main, xlims(1) + 0.05*diff(xlims), ylims(2) - 0.1*diff(ylims), text_str, ...
    'FontSize', 12, 'EdgeColor', 'black', 'BackgroundColor', 'white');
% --- SAVE FIGURE ---
print(fig3, 'Adsorption/LineAverageModel/Figures/figure_spr_shift_composite.png', '-dpng', '-r300');
print(fig3, 'Adsorption/LineAverageModel/Figures/EPS/figure_spr_shift_composite.eps', '-depsc');
%% --- STEP 5: VIDEO CREATION ---
% =========================================================================
if generate_video_frames
    fprintf('\n--- Iniciando a criação do vídeo a partir dos frames salvos ---\n');
    tic;
    % --- THE FIX IS HERE (Part 1): Change the video profile and filename ---
    video_filename = 'Adsorption/LineAverageModel/Videos/spri_simulation_final.avi';
    outputVideo = VideoWriter(video_filename, 'Motion JPEG AVI');
    
    outputVideo.FrameRate = 30;
    outputVideo.Quality = 95; % Quality for AVI is 0-100
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
fig3b = figure('Name', 'Proper Sensorgram: Resonance Angle vs. Time', 'Position', [300, 300, 1400, 700]);
lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);
for i = 1:length(lines_to_plot)
    subplot(1, length(lines_to_plot), i);
    line_idx = lines_to_plot(i);
    
    % This is the proper, physically-correct sensorgram
    plot(t_exp, theta_spr_vs_time(:, line_idx), 'r-', 'LineWidth', 2, 'DisplayName', 'Sensorgram (Resonance Angle)');
    
    grid on;
    xlabel('Time (s)');
    ylabel('Resonance Angle (degrees)');
    title(sprintf('Sensorgram for Line %d', line_idx));
    
    % For validation, plot the original RU data on a separate y-axis
    yyaxis right
    plot(t_exp, s_obs_ru(:, line_idx), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Original Simulation (RU)');
    ylabel('Response Units (RU)');
    
    legend('Location', 'best');
    xlim([0, t_exp(end)]);
end
%sgtitle('Final Validation: The Physically Correct Sensorgram (Angle vs. Time)', 'FontSize', 16);
% --- SAVE FIGURE ---
print(fig3b, 'Adsorption/LineAverageModel/Figures/figure_angle_vs_ru.png', '-dpng', '-r300');
print(fig3b, 'Adsorption/LineAverageModel/Figures/EPS/figure_angle_vs_ru.eps', '-depsc');
%% 
global_sensorgram = sum(formula_response_vs_time, 2);
% --- Step 2: Create the Figure and Subplots ---
fig_formula = figure('Name', 'Final Sensorgram Results from Formula', 'Position', [200, 200, 1400, 600]);
% --- Plot 1: All Individual Line Sensorgrams ---
subplot(1, 2, 1);
plot(t_exp, formula_response_vs_time, 'LineWidth', 1.5);
grid on;
title('Individual Line Sensorgrams');
xlabel('Time (s)');
ylabel('Change from Baseline (\Delta{N}_s^{eff})');
xlim([0, t_exp(end)]);
% Optional: Add a legend if you have a small number of lines
if ads_y_dim <= 10
    legend(arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), 'Location', 'best');
end
% --- Plot 2: Global (Summed) Sensorgram ---
subplot(1, 2, 2);
plot(t_exp, global_sensorgram, 'r-', 'LineWidth', 2);
grid on;
title('Global (Summed) Sensorgram');
xlabel('Time (s)');
ylabel('Total Change (\Sigma \Delta{N}_s^{eff})');
xlim([0, t_exp(end)]);
sgtitle('Final Sensorgrams Calculated from Analytical Formula', 'FontSize', 16, 'FontWeight', 'bold');
% --- SAVE FIGURE ---
print(fig_formula, 'Adsorption/LineAverageModel/Figures/figure_final_sensorgrams_formula.png', '-dpng', '-r300');
print(fig_formula, 'Adsorption/LineAverageModel/Figures/EPS/figure_final_sensorgrams_formula.eps', '-depsc');

absolute_neff_vs_time = formula_response_vs_time + n_bulk;
% --- Step 2: Calculate the Global (Summed) Absolute N_s_eff ---
global_absolute_neff = n_bulk+global_sensorgram;
% --- Step 3: Create the Figure and Subplots ---
fig_abs_ri = figure('Name', 'Absolute Effective RI Sensorgrams', 'Position', [200, 200, 1400, 600]);
% --- Plot 1: All Individual Line Sensorgrams ---
subplot(1, 2, 1);
plot(t_exp, absolute_neff_vs_time, 'LineWidth', 1.5);
grid on;
title('Individual Line Sensorgrams (Absolute N_s^{eff})');
xlabel('Time (s)');
ylabel('Absolute Effective RI ({N}_s^{eff})');
xlim([0, t_exp(end)]);
% Optional: Add a legend if you have a small number of lines
if ads_y_dim <= 10
    legend(arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), 'Location', 'best');
end
% --- Plot 2: Global (Summed) Sensorgram ---
subplot(1, 2, 2);
plot(t_exp, global_absolute_neff, 'r-', 'LineWidth', 2);
grid on;
title('Global (Summed) Sensorgram');
xlabel('Time (s)');
ylabel('Absolute Effective RI (\Sigma {N}_s^{eff})');
xlim([0, t_exp(end)]);
sgtitle('Final Sensorgrams Plotted as Absolute N_s^{eff}', 'FontSize', 16, 'FontWeight', 'bold');
% --- SAVE FIGURE ---
print(fig_abs_ri, 'Adsorption/LineAverageModel/Figures/figure_absolute_ri_sensorgrams.png', '-dpng', '-r300');
print(fig_abs_ri, 'Adsorption/LineAverageModel/Figures/EPS/figure_absolute_ri_sensorgrams.eps', '-depsc');
% =========================================================================
%% --- STEP 7: INVERSE PROCESS ---
% =========================================================================
% --- Step 7.1: Analyze saved frames to reconstruct the sensorgram ---
theta_extracted_from_frames = analyze_spr_frames_to_get_sensorgram(video_frames_folder, angle_range);
% --- Step 7.2: Convert extracted angles back to Refractive Index ---
fprintf('Converting extracted angles back to Refractive Index via interpolation...\n');
n2_reconstructed_vs_time = zeros(size(theta_extracted_from_frames));
parfor j_idx = 1:ads_y_dim
    [unique_thetas, unique_indices] = unique(theta_spr_vs_time(:, j_idx));
    unique_n2s = n2_vs_time(unique_indices, j_idx);
    n2_reconstructed_vs_time(:, j_idx) = interp1(unique_thetas, unique_n2s, theta_extracted_from_frames(:, j_idx), 'linear', 'extrap');
end
% --- Step 7.3: Convert reconstructed RI back to RU ---
ru_reconstructed = (n2_reconstructed_vs_time - n_bulk) / RU_TO_RIU;
% =====================================================================
% --- FINAL VALIDATION PLOTS WITH QUANTITATIVE ERROR IN LEGENDS ---
% =====================================================================
lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);

% --- Plot A: Validation in Resonance Angle units ---
fig4a = figure('Name', 'Final Validation: Original vs. Extracted Angle', 'Position', [300, 300, 1800, 500]);
ax_handles_a = gobjects(1, length(lines_to_plot)); % Preallocate handles array
for i = 1:length(lines_to_plot)
    ax_handles_a(i) = subplot(1, length(lines_to_plot), i);
    line_idx = lines_to_plot(i);
    original_angle = theta_spr_vs_time(:, line_idx);
    extracted_angle = theta_extracted_from_frames(:, line_idx);
    
    [offset, error_abs] = calculate_validation_error(original_angle, extracted_angle);
    extracted_legend_text = sprintf('Extracted (Offset=%.1e, Error=%.1e)', offset, error_abs);
    plot(t_exp, original_angle - original_angle(1), 'b-', 'LineWidth', 4, 'DisplayName', 'Original');
    hold on;
    plot(t_exp, extracted_angle - extracted_angle(1), 'r--', 'LineWidth', 2, 'DisplayName', extracted_legend_text);
    
    grid on; xlabel('Time (s)'); ylabel('Change in Resonance Angle (degrees)');
    title(sprintf('Angle Validation for Line %d', line_idx));
    legend('Location', 'best'); xlim([0, t_exp(end)]);
end
%sgtitle('Final Validation in Angle Space', 'FontSize', 16);
% --- Save full figure ---
print(fig4a, 'Adsorption/LineAverageModel/Figures/figure_4a_validation_angle.png', '-dpng', '-r300');
print(fig4a, 'Adsorption/LineAverageModel/Figures/EPS/figure_4a_validation_angle.eps', '-depsc');
% --- Save individual subplots ---
fprintf('Saving individual angle validation subplots as EPS files...\n');
base_path = 'Adsorption/LineAverageModel/Figures/EPS/';
for k = 1:length(ax_handles_a)
    line_idx = lines_to_plot(k);
    filename = sprintf('%svalidation_angle_line_%d.eps', base_path, line_idx);
    save_subplot_as_eps(ax_handles_a(k), filename);
end
fprintf('Finished saving angle validation subplots.\n');


% --- Plot B: Validation in Refractive Index units ---
fig4b = figure('Name', 'Final Validation: Original vs. Reconstructed RI', 'Position', [300, 300, 1800, 500]);
ax_handles_b = gobjects(1, length(lines_to_plot)); % Preallocate handles array
for i = 1:length(lines_to_plot)
    ax_handles_b(i) = subplot(1, length(lines_to_plot), i);
    line_idx = lines_to_plot(i);
    original_ri = n2_vs_time(:, line_idx);
    reconstructed_ri = n2_reconstructed_vs_time(:, line_idx);
    
    [offset, error_abs] = calculate_validation_error(original_ri, reconstructed_ri);
    extracted_legend_text = sprintf('Reconstructed (Offset=%.1e, Error=%.1e)', offset, error_abs);
    plot(t_exp, original_ri - original_ri(1), 'b-', 'LineWidth', 4, 'DisplayName', 'Original');
    hold on;
    plot(t_exp, reconstructed_ri - reconstructed_ri(1), 'r--', 'LineWidth', 2, 'DisplayName', extracted_legend_text);
    
    grid on; xlabel('Time (s)'); ylabel('Change in Refractive Index (dRIU)');
    title(sprintf('RI Validation for Line %d', line_idx));
    legend('Location', 'best'); xlim([0, t_exp(end)]);
end
%sgtitle('Final Validation in Refractive Index Space', 'FontSize', 16);
% --- Save full figure ---
print(fig4b, 'Adsorption/LineAverageModel/Figures/figure_4b_validation_ri.png', '-dpng', '-r300');
print(fig4b, 'Adsorption/LineAverageModel/Figures/EPS/figure_4b_validation_ri.eps', '-depsc');
% --- Save individual subplots ---
fprintf('Saving individual RI validation subplots as EPS files...\n');
for k = 1:length(ax_handles_b)
    line_idx = lines_to_plot(k);
    filename = sprintf('%svalidation_ri_line_%d.eps', base_path, line_idx);
    save_subplot_as_eps(ax_handles_b(k), filename);
end
fprintf('Finished saving RI validation subplots.\n');


% --- Plot C: Validation in Response Units ---
fig4c = figure('Name', 'Final Round-Trip Validation: RU Original vs. Reconstructed', 'Position', [300, 300, 1800, 500]);
ax_handles_c = gobjects(1, length(lines_to_plot)); % Preallocate handles array
for i = 1:length(lines_to_plot)
    ax_handles_c(i) = subplot(1, length(lines_to_plot), i);
    line_idx = lines_to_plot(i);
    original_ru = s_obs_ru(:, line_idx);
    reconstructed_ru_line = ru_reconstructed(:, line_idx);
    
    [offset, error_abs] = calculate_validation_error(original_ru, reconstructed_ru_line);
    extracted_legend_text = sprintf('Reconstructed (Offset=%.1e, Error=%.1e)', offset, error_abs);
    plot(t_exp, original_ru - original_ru(1), 'b-', 'LineWidth', 4, 'DisplayName', 'Original');
    hold on;
    plot(t_exp, reconstructed_ru_line - reconstructed_ru_line(1), 'r--', 'LineWidth', 2, 'DisplayName', extracted_legend_text);
    
    grid on; xlabel('Time (s)'); ylabel('Response Units (RU)');
    title(sprintf('RU Validation for Line %d', line_idx));
    legend('Location', 'best'); xlim([0, t_exp(end)]);
end
%sgtitle('Final Round-Trip Validation in Response Units', 'FontSize', 16);
% --- Save full figure ---
print(fig4c, 'Adsorption/LineAverageModel/Figures/figure_4c_validation_ru.png', '-dpng', '-r300');
print(fig4c, 'Adsorption/LineAverageModel/Figures/EPS/figure_4c_validation_ru.eps', '-depsc');
% --- Save individual subplots ---
fprintf('Saving individual RU validation subplots as EPS files...\n');
for k = 1:length(ax_handles_c)
    line_idx = lines_to_plot(k);
    filename = sprintf('%svalidation_ru_line_%d.eps', base_path, line_idx);
    save_subplot_as_eps(ax_handles_c(k), filename);
end
fprintf('Finished saving RU validation subplots.\n');
% =========================================================================
% --- FINAL PUBLICATION PLOT: ANNOTATED SENSORGRAM FIT ---
% =========================================================================
fprintf('\nGenerating annotated experiment plot for publication...\n');
% --- Select data to plot (e.g., the first experiment and the middle line) ---
exp_to_plot = 1;
line_to_plot = round(ads_y_dim/2);
% Get the relevant data from your previous calculations
setting_to_plot = exp_settings(exp_to_plot);
t_exp = exp_data{exp_to_plot}.time;
% Get the noisy experimental data for the chosen line
s_obs_noisy = exp_data{exp_to_plot}.signals(:, line_to_plot);
% Reconstruct the final smooth fit from the optimized parameters
% We need to run the model one last time with the final 'opt_params_1D'
[~, s_final_fit_all_lines] = run_single_experiment_1D_model(opt_params_1D, setting_to_plot, model_config);
s_final_fit = s_final_fit_all_lines(:, line_to_plot);
% --- Create the Plot ---
fig_annotated = figure('Name', 'Annotated Sensorgram Fit', 'Position', [100, 100, 1000, 600]);
hold on;
% Plot the noisy data and the final smooth fit
plot(t_exp, s_obs_noisy, '.', 'Color', [0.6 0.6 1], 'DisplayName', 'Noisy Experimental Data'); % Light blue dots for data
plot(t_exp, s_final_fit, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Final Model Fit'); % Bold red line for fit
% --- Add vertical lines and text annotations for each phase ---
y_lims = ylim; % Get current y-axis limits
% Set text position to be 95% of the axis height (inside the plot)
text_y_pos = y_lims(1) + 0.95 * (y_lims(2) - y_lims(1)); 
% Get all time breaks and concentrations for annotation
all_t_breaks = [0, setting_to_plot.pulse_times, t_total];
all_concs = [setting_to_plot.pulse_concs, c_diss];
for i = 1:length(all_concs)
    t_start = all_t_breaks(i);
    t_end = all_t_breaks(i+1);
    
    % Draw vertical line at the start of the new phase
    if i > 1
        line([t_start, t_start], y_lims, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    % Add text annotation in the middle of the phase
    text_x_pos = (t_start + t_end) / 2;
    text_str = sprintf('C = %.1e M', all_concs(i));
    
    % Add a background to the text for better legibility
    text(text_x_pos, text_y_pos, text_str, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [1 1 1 0.7], ... % White, 70% opaque background
        'Margin', 3); % Padding around text
end
% --- Finalize Plot ---
hold off;
box on;
grid on;
axis tight; % Ensure plot fits data snugly
% title(sprintf('Final Model Fit to Noisy Data for Line %d', line_to_plot));
xlabel('Time (s)');
ylabel('Response (RU)');
legend('Location', 'southeast');
% --- SAVE FIGURE ---
print(fig_annotated, 'Adsorption/LineAverageModel/Figures/figure_6_annotated_fit.png', '-dpng', '-r300');
print(fig_annotated, 'Adsorption/LineAverageModel/Figures/EPS/figure_6_annotated_fit.eps', '-depsc');



