%% SCRIPT FOR INVERSE PROCESSING OF EXPERIMENTAL SPRI VIDEO
% MODIFIED TO USE DIRECT ANALYTICAL INVERSE FORMULA

clearvars; close all; clc;

% =========================================================================
%% --- STEP 1: USER CONFIGURATION & PHYSICAL CONSTANTS ---
% =========================================================================
% (This section remains unchanged)
video_filename = 'Adsorption/LineAverageModel/Videos/spri_simulation_final.avi';
wavelength = 670; d1 = 50; d2 = 1000;
n0 = sqrt(2.3104);
n1 = sqrt(-14.379 + 1.0084j);
n_bulk = sqrt(1.7876);
RU_TO_RIU = 0.001 / 1000;
angle_range = linspace(65, 80, 1280);

fprintf('Configuration loaded for video: %s\n', video_filename);

% =========================================================================
%% --- STEP 2: READ VIDEO AND EXTRACT FRAMES ---
% =========================================================================
% (This section remains unchanged)
fprintf('\n--- Step 2: Reading video file and extracting frames... ---\n');
tic;
try
    videoObj = VideoReader(video_filename);
    fprintf('✅ VideoReader successfully created the object.\n');
catch ME
    fprintf(2, '\n--- DETAILED MATLAB ERROR INFORMATION ---\n');
    fprintf(2, 'Error Message:     %s\n', ME.message);
    fprintf(2, 'Error Identifier:  %s\n', ME.identifier);
    if ~isempty(ME.stack), fprintf(2, 'Error occurred in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line); end
    fprintf(2, '--- END OF ERROR --- \n\n');
    error('Execution stopped because the video file could not be read.');
end
numFrames = videoObj.NumFrames;
vidHeight = videoObj.Height;
vidWidth = videoObj.Width;
frameRate = videoObj.FrameRate;
if vidWidth ~= length(angle_range), error('Video width (%d) does not match angle_range length (%d).', vidWidth, length(angle_range)); end
fprintf('Video details: %d frames, %d lines (height), %d angle points (width).\n', numFrames, vidHeight, vidWidth);
all_frames = zeros(vidHeight, vidWidth, numFrames, 'uint8');
frame_count = 0;
while hasFrame(videoObj)
    frame_count = frame_count + 1;
    frame_rgb = readFrame(videoObj);
    if size(frame_rgb, 3) == 3, all_frames(:, :, frame_count) = rgb2gray(frame_rgb); else, all_frames(:, :, frame_count) = frame_rgb; end
    if mod(frame_count, 100) == 0, fprintf('Processed frame %d of %d...\n', frame_count, numFrames); end
end
fprintf('Successfully read %d frames.\n', frame_count);
toc;

% =========================================================================
%% --- STEP 3 (REVISED): ROBUST ANGLE EXTRACTION WITH PARABOLIC FIT ---
% =========================================================================
fprintf('\n--- Step 3 (Revised): Extracting angles using robust parabolic fit... ---\n');
tic;

theta_vs_time = zeros(numFrames, vidHeight);

% Define the size of the window for the parabolic fit (+/- pixels from the rough minimum)
fit_window_half_size = 130; 
for t_idx = 1:numFrames
    current_frame = double(all_frames(:, :, t_idx)); 

    for line_idx = 1:vidHeight
        spr_curve = current_frame(line_idx, :);

        % --- REFINEMENT: Spatially smooth the curve to reduce pixel noise before fitting ---
        % Find the approximate minimum on the SMOOTHED curve
        [~, rough_min_idx] = min(spr_curve);

        % Define the window on the SMOOTHED curve for fitting
        start_idx = max(1, rough_min_idx - fit_window_half_size);
        end_idx   = min(length(spr_curve), rough_min_idx + fit_window_half_size);

        if (end_idx - start_idx) < 2 % Need at least 3 points for a parabolic fit
            theta_vs_time(t_idx, line_idx) = angle_range(rough_min_idx);
            continue;
        end

        x_fit = start_idx:end_idx;
        y_fit = spr_curve(x_fit); % Use the smoothed data for the fit

        % Fit a 2nd-degree polynomial (parabola). This is very fast.
        p = polyfit(x_fit, y_fit, 2);

        % Find the analytical minimum of the parabola: x = -b / (2a)
        subpixel_min_idx = -p(2) / (2 * p(1));
        
        % Convert the sub-pixel index back to a real angle
        theta_vs_time(t_idx, line_idx) = interp1(1:length(angle_range), angle_range, subpixel_min_idx);
    end
end

fprintf('Refined angle extraction complete.\n');
toc;
%% <<< THIS SECTION REPLACES THE OLD STEP 4 AND 5 >>>
% =========================================================================
%% --- STEP 4 (REVISED): DIRECT INVERSE CALCULATION FROM ANGLE ---
% =========================================================================
fprintf('\n--- Step 4 (Revised): Directly calculating RI from angles using analytical formula... ---\n');
tic;

% Define constants from the formula using your script's variables
% n2 in the formula is the prism's refractive index (n0 in your script)
n_prism = n0; 
% ε_mr in the formula is the real part of the metal's dielectric constant
epsilon_metal_complex = n1^2;
epsilon_mr = real(epsilon_metal_complex);

% Convert the measured angles from degrees to radians for sin()
theta_rad = deg2rad(theta_vs_time);

% --- Apply the formula ---
% Calculate the term (n_prism * sin(θ_res))^2
term_sin_sq = (n_prism * sin(theta_rad)).^2;

% Calculate the numerator and denominator of the main fraction
numerator = epsilon_mr * term_sin_sq;
denominator = epsilon_mr - term_sin_sq;

% Calculate n3^2, which we call n_eff^2
n_eff_sq = numerator ./ denominator;

% n3 is the square root, which is our effective refractive index (n_eff)
n_eff = sqrt(n_eff_sq);

% --- Convert the calculated effective refractiv
% e index to Response Units ---
% The signal in RU is proportional to the change from the bulk refractive index.
% --- This is the end of your revised Step 4 ---
ru_reconstructed = (n_eff - n_bulk) / RU_TO_RIU;

fprintf('Direct calculation and conversion to RU complete.\n');
toc;

% =========================================================================
%% --- STEP 5: PLOT RECONSTRUCTED SENSORGRAMS ---
% =========================================================================
fprintf('\n--- Step 5: Plotting results... ---\n');
time_axis = (0:numFrames-1) / frameRate;

figure('Name', 'Reconstructed Sensorgrams (Direct Formula, Corrected)', 'Position', [100, 100, 900, 600]);
plot(time_axis, ru_reconstructed, 'LineWidth', 1.5);

grid on; box on;
title('Reconstructed Sensorgrams (Baseline Corrected)');
xlabel('Time (s)');
ylabel('Response (RU)');
xlim([0, time_axis(end)]);
legend_labels = arrayfun(@(j) sprintf('Line %d', j), 1:vidHeight, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northwest');
fprintf('Plotting complete.\n');

% =========================================================================
%% --- STEP 6: DEFINE EXPERIMENTAL CONDITIONS ---
% =========================================================================
% ⚠️ ACTION REQUIRED: You must update these values to match the
% conditions used to collect your experimental data.

fprintf('Defining experimental conditions...\n');
model_config.D_coeff = 6e-5;
model_config.ru_to_m = 1e-10;
gridN_x = 22; gridN_y = 5; gridN_z = 3;ads_layer = 1;
ads_x_range = [5,15]; ads_y_range = [1,5];
% Get number of lines in adsorption region
ads_y_dim = ads_y_range(2) - ads_y_range(1) + 1;
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
% You may need to add other model_config fields if your functions require them
% model_config.gridN_x = 22; etc...


% =========================================================================
%% --- STEP 7: SETUP OPTIMIZER (Initial Guess & Bounds) ---
% =========================================================================
% This section is adapted from your original script to set up the optimizer.
fprintf('Setting up optimizer with initial guesses and bounds...\n');

% Initial guess based on typical homogeneous parameters
homog_params = [9.4e3, 0.0078, 2960]; % [kon, koff, smax_total]
p0_base_kon = log10(homog_params(1)) * ones(ads_y_dim, 1);
p0_base_koff = log10(homog_params(2)) * ones(ads_y_dim, 1);
p0_base_smax = log10(homog_params(3) / ads_y_dim) * ones(ads_y_dim, 1);
p0_base = [p0_base_kon; p0_base_koff; p0_base_smax];

% Define parameter bounds in log10 space to keep the search stable
lb_kon = log10(1e2); ub_kon = log10(1e7);
lb_koff = log10(1e-6); ub_koff = log10(1e-1);
lb_smax = log10(1); ub_smax = log10(5000);
lb = [repmat(lb_kon, ads_y_dim, 1); repmat(lb_koff, ads_y_dim, 1); repmat(lb_smax, ads_y_dim, 1)];
ub = [repmat(ub_kon, ads_y_dim, 1); repmat(ub_koff, ads_y_dim, 1); repmat(ub_smax, ads_y_dim, 1)];

% Create the final initial guess (p0) by adding a small random perturbation
rng('default');
p0 = max(lb, p0_base + 0.1 * randn(size(p0_base)));
p0 = min(ub, p0);
[kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
    create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
kon_ads_2D = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
koff_ads_2D = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
smax_ads_2D = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
% Average parameters along the x-axis (flow direction) to get line parameters
p_true_kon_1D = mean(kon_ads_2D, 1)';   % [ads_y_dim x 1]
p_true_koff_1D = mean(koff_ads_2D, 1)'; % [ads_y_dim x 1]
p_true_smax_1D = mean(smax_ads_2D, 1)'; % [ads_y_dim x 1]
p_true_1D = [p_true_kon_1D; p_true_koff_1D; p_true_smax_1D];
N_params_1D = length(p_true_1D);
exp_data = cell(1, 1);
data_struct.time = time_axis;
data_struct.signals = ru_reconstructed; % Your experimental data
exp_data{1} = data_struct;
optim_plot_fun = @(log_p, optim_v, state) optimPlotter_1D(...
    log_p, optim_v, state, ...
    p_true_1D, exp_data, ads_y_dim);
% Define optimization options
optim_opts = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter', ...
    'MaxIterations', 10, ... % Increased iterations for real data
    'UseParallel', true, ...
    'FunctionTolerance', 1e-8, ...
    'StepTolerance', 1e-8, ...
    'OutputFcn', optim_plot_fun);
    

% =========================================================================
%% --- STEP 8: RUN PARAMETER IDENTIFICATION ---
% =========================================================================
% Prepare the data in the format expected by the residual function


% Define the function handle for the optimizer
% This tells lsqnonlin what function to minimize
residual_fun = @(log_params) compute_residuals_1D_model(log_params, exp_settings, exp_data, model_config);

fprintf('\nStarting parameter identification...\n');
tic;

% Run the optimization
[opt_log_params, resnorm] = lsqnonlin(residual_fun, p0, lb, ub, optim_opts);

toc;
fprintf('Parameter identification complete.\n');

% Convert final parameters from log10 scale to linear scale
opt_params_1D = 10.^opt_log_params;


% =========================================================================
%% --- STEP 9: DISPLAY AND PLOT FINAL RESULTS ---
% =========================================================================
% --- Display Final Parameters ---
fprintf('\n--- Optimal Kinetic Parameters ---\n');
opt_kon = opt_params_1D(1:ads_y_dim);
opt_koff = opt_params_1D(ads_y_dim+1 : 2*ads_y_dim);
opt_smax = opt_params_1D(2*ads_y_dim+1 : 3*ads_y_dim);

param_table = table((1:ads_y_dim)', opt_kon, opt_koff, opt_smax, ...
    'VariableNames', {'Line', 'k_on', 'k_off', 's_max'});
disp(param_table);

% --- Plot Final Fit vs. Experimental Data ---
fprintf('\nGenerating final plot of model fit vs. experimental data...\n');

% Simulate the model one last time with the optimal parameters
[~, s_final_fit] = run_single_experiment_1D_model(opt_params_1D, exp_settings, model_config);

figure('Name', 'Model Fit vs. Experimental Data', 'Position', [100, 100, 1200, 700]);
num_lines = ads_y_dim;
for i = 1:num_lines
    subplot(ceil(num_lines/2), 2, i);
    plot(time_axis, ru_reconstructed(:, i), 'b.', 'MarkerSize', 4, 'DisplayName', 'Experimental Data');
    hold on;
    plot(time_axis, s_final_fit(:, i), 'r-', 'LineWidth', 2, 'DisplayName', 'Model Fit');
    hold off;
    grid on;
    title(sprintf('Line %d', i));
    xlabel('Time (s)');
    ylabel('Response (RU)');
    legend('Location', 'best');
end
sgtitle('Final Model Fit Compared to Experimental Data', 'FontSize', 14, 'FontWeight', 'bold');