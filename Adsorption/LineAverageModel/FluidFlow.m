%% Fluid Dynamics and Adsorption Visualization Script
% This script runs a full 3D forward simulation to visualize the analyte
% concentration in the flow cell and the corresponding adsorption on the
% sensor surface over time.
%
% Assumes the following functions are in the MATLAB path:
% - simulate_3d_flow_model_with_pulses
% - create_ground_truth_heterogeneity
% - create_velocity_profile
% - generate_experiments
%
clearvars; close all; clc;

fprintf('--- Setting up Dynamic Visualization ---\n');

% --- Create the Videos folder if it doesn't exist ---
video_folder_path = 'Adsorption/LineAverageModel/Videos';
if ~exist(video_folder_path, 'dir')
    mkdir(video_folder_path);
end

%% --- STEP 1: Define Physical and Experimental Parameters ---
% These parameters should match your main script to ensure consistency.
% ========================================================================
% --- Grid and Adsorption Region ---
gridN_x = 22; gridN_y = 5; gridN_z = 15;
ads_layer = 1;
ads_x_range = [5, 15]; ads_y_range = [1, 5];

% --- Physical Constants ---
D_coeff = 6e-5;
ru_to_m = 1e-10;
grid_size_x = 11.0; % mm
grid_size_z = 0.3;  % mm
dx = grid_size_x / gridN_x;
dz = grid_size_z / gridN_z;

% --- Experiment Parameters (single experiment for visualization) ---
n_exp = 1;
base_max_velocity = 8.3;
T1 = 2200; T2 = 2 * T1; T3 = 3 * T1; t_total = 4 * T1;
c_diss = 0; c1 = 3.3e-4; c2 = 0.21e-4;
exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);
setting = exp_settings(1); % Use the first (and only) experiment setting

% --- Create Ground Truth Heterogeneous Surface ---
[kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
    create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);

% --- Create Velocity Profile ---
velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);

fprintf('Parameters loaded and ground truth created.\n');

%% --- STEP 2: Run Full 3D Simulation to Get Time-Dependent Data ---
% This is the core computational step. We capture the full 3D matrices
% for concentration (C) and surface coverage (s) at each time point.
% ========================================================================
fprintf('\nRunning full 3D simulation... This may take a moment.\n');
tic;
s0_grid = zeros(gridN_x, gridN_y, gridN_z);
t_breaks = [0, setting.pulse_times, setting.t_total];
concentrations = [setting.pulse_concs, setting.c_diss];

% *** CRUCIAL: Capture all outputs, especially C_grid_t and s_grid_t ***
[t, C_grid_t, s_grid_t] = simulate_3d_flow_model_with_pulses_fluid(...
    gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
    velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid, dx, dz);

fprintf('Simulation complete. Generated %d time points.\n', length(t));
toc;

%% --- STEP 3: Create Dynamics Animation ---
% We will now loop through each time step and create a video frame showing
% the state of the system.
% ========================================================================
fprintf('\nPreparing to generate animation...\n');

% --- Video Writer Setup ---
video_filename = fullfile(video_folder_path, 'flow_and_adsorption_dynamics.avi');
outputVideo = VideoWriter(video_filename, 'Motion JPEG AVI');
outputVideo.FrameRate = 10; % Adjust for desired speed
outputVideo.Quality = 100;
open(outputVideo);

% --- Figure and Axes Setup ---
fig = figure('Name', 'Flow and Adsorption Dynamics', 'Position', [50, 50, 1400, 600]);
tlo = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Left Plot: Fluid Concentration (Side View) ---
ax1 = nexttile;
% Take a slice through the middle of the y-dimension to get an x-z view
C_slice_xz = squeeze(C_grid_t(:, round(gridN_y/2), :, 1));
h_conc = imagesc(ax1, C_slice_xz');
ax1.YDir = 'normal'; % Place z=0 at the bottom
title(ax1, 'Analyte Concentration in Flow Cell (Side View)');
xlabel(ax1, 'Position along Flow (x-index)');
ylabel(ax1, 'Position in Height (z-index)');
colorbar(ax1);
caxis(ax1, [0, max(concentrations)]); % Fix color axis for consistency

% --- Right Plot: Surface Adsorption (Top View) ---
ax2 = nexttile;
% Take a slice at the sensor layer (z = ads_layer) to get an x-y view
s_slice_xy = squeeze(s_grid_t(:, :, ads_layer, 1));
h_surf = imagesc(ax2, s_slice_xy');
ax2.YDir = 'normal';
title(ax2, 'Bound Analyte on Sensor Surface (Top View)');
xlabel(ax2, 'Position along Flow (x-index)');
ylabel(ax2, 'Line Index (y-index)');
colorbar(ax2);
caxis(ax2, [0, max(s_grid_t(:)) + eps]); % Fix color axis

% --- Animation Loop ---
fprintf('Generating %d frames for the video...\n', length(t));
tic;
for i = 1:length(t)
    % --- Extract data for the current time step ---
    current_C_grid = C_grid_t(:, :, :, i);
    current_s_grid = s_grid_t(:, :, :, i);

    % --- Prepare the 2D slices for plotting ---
    C_slice_xz = squeeze(current_C_grid(:, round(gridN_y/2), :));
    s_slice_xy = squeeze(current_s_grid(:, :, ads_layer));

    % --- Update the data in the existing plots (more efficient) ---
    set(h_conc, 'CData', C_slice_xz');
    set(h_surf, 'CData', s_slice_xy');

    % --- Update the main title with the current time ---
    sgtitle(tlo, sprintf('System Dynamics at Time = %.1f s', t(i)), 'FontWeight', 'bold');

    % --- Capture the figure as a video frame ---
    frame = getframe(fig);
    writeVideo(outputVideo, frame);

    % --- Display progress ---
    if mod(i, 50) == 0
        fprintf('  ... processed frame %d of %d\n', i, length(t));
    end
end
toc;

% --- Finalize and Close ---
close(outputVideo);
close(fig);
fprintf('\nAnimation saved successfully to:\n%s\n', video_filename);