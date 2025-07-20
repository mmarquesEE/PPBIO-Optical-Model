%% Main Script to Generate a High-Quality SPR Visualization Video
% This script uses a user-provided 4-layer Fresnel function to link the
% simulated sensorgram to the underlying optical phenomena.
%
clearvars; close all; clc;

fprintf('--- Setting up SPR Physics Visualization Video ---\n');

%% --- STEP 1: Define Physical and Experimental Parameters ---
% ========================================================================
gridN_x = 22; gridN_y = 5; gridN_z = 15;
ads_layer = 1;
ads_x_range = [5, 15]; ads_y_range = [1, 5];
D_coeff = 6e-5;
ru_to_m = 1e-10;
grid_size_x = 11.0; 
grid_size_z = 0.3;  
dx = grid_size_x / gridN_x;
dz = grid_size_z / gridN_z;
n_exp = 1;
base_max_velocity = 8.3e-2;
T1 = 2200; T2 = 2 * T1; T3 = 3 * T1; t_total = 4 * T1;
c_diss = 0; c1 = 3.3e-6; c2 = 0.21e-6;
exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);
setting = exp_settings(1); 
[kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
    create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);
fprintf('Parameters loaded and ground truth created.\n');

%% --- STEP 2: Run Full 3D Simulation ---
% ========================================================================
fprintf('\nRunning full 3D simulation... This may take a moment.\n');
tic;
s0_grid = zeros(gridN_x, gridN_y, gridN_z);
t_breaks = [0, setting.pulse_times, setting.t_total];
concentrations = [setting.pulse_concs, setting.c_diss];

[t, ~, s_grid_t] = simulate_3d_flow_model_with_pulses_fluid(...
    gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
    velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid, dx, dz);
fprintf('Simulation complete. Generated %d time points.\n', length(t));
toc;

%% --- STEP 3: Pre-process Data and Define Optical Parameters ---
% ========================================================================
fprintf('Pre-processing data for animation...\n');

% Calculate Line-Level Sensorgrams (this is our "RU" data)
s_ads_region_t = squeeze(s_grid_t(ads_x_range(1):ads_x_range(2), :, ads_layer, :));
line_grams_ru = squeeze(sum(s_ads_region_t, 1)); % Dims: [ny, time]

% --- Define Optical and Physical Constants ---
wavelength = 670; % nm
d1 = 50;          % Gold film thickness (nm)
d2 = 1000;        % Analyte layer thickness (nm), kept large as in your example
n0 = sqrt(2.3104);         % Layer 0: Optical substrate (Prism)
n1 = sqrt(-14.379 + 1.0084j); % Layer 1: Gold film (complex RI)
n_bulk = sqrt(1.7876);     % Layer 3: Flow cell solution (baseline buffer)
angle_range = linspace(65, 80, 1080); % Angular range for SPR curve
RU_TO_RIU = 0.001 / 1000; % Conversion factor

%% --- STEP 4: Call the SPR Animation Function ---
% ========================================================================
create_spr_animation(t, line_grams_ru, t_breaks, ...
                     wavelength, d1, d2, n0, n1, n_bulk, angle_range, RU_TO_RIU);

fprintf('\nSPR physics video animation saved successfully.\n');


%% ========================================================================
%  --- LOCAL FUNCTIONS ---
%  ========================================================================
function create_spr_animation(t, line_grams_ru, t_breaks, wavelength, d1, d2, n0, n1, n_bulk, angle_range, RU_TO_RIU)
    
    filename = 'spr_physics_animation.mp4';
    outputVideo = VideoWriter(filename, 'MPEG-4');
    outputVideo.Quality = 100;
    outputVideo.FrameRate = 30;
    open(outputVideo);

    background_color = [217/255, 217/255, 217/255];
    fig = figure('Position', [50 50 1920 800], 'Color', background_color);
    
    tlo = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    num_lines = size(line_grams_ru, 1);
    line_colors = lines(num_lines);

    % --- Setup Panel 1: Line Sensorgrams ---
    ax1 = nexttile;
    hold(ax1, 'on');
    animated_lines = gobjects(num_lines, 1);
    for k = 1:num_lines, animated_lines(k) = animatedline(ax1, 'Color', line_colors(k,:), 'LineWidth', 2); end
    hold(ax1, 'off'); grid(ax1, 'on');
    title(ax1, 'Line Sensorgrams', 'FontSize', 14);
    xlabel(ax1, 'Time (s)'); ylabel(ax1, 'Response (RU)');
    xlim(ax1, [0, t(end)]);
    legend(ax1, arrayfun(@(k) sprintf('Line %d', k), 1:num_lines, 'UniformOutput', false), 'Location', 'northwest');
    time_marker1 = xline(ax1, 0, 'r--', 'LineWidth', 1.5);

    % --- Setup Panel 2: SPR Reflectance Curves ---
    ax2 = nexttile;
    hold(ax2, 'on');
    reflectance_plots = gobjects(num_lines, 1);
    for k = 1:num_lines, reflectance_plots(k) = plot(ax2, angle_range, zeros(size(angle_range)), 'Color', line_colors(k,:), 'LineWidth', 2); end
    hold(ax2, 'off'); grid(ax2, 'on');
    title(ax2, 'SPR Reflectance Curves', 'FontSize', 14);
    xlabel(ax2, 'Incident Angle (degrees)'); ylabel(ax2, 'Reflectivity');
    ylim(ax2, [0, 1]); xlim(ax2, [angle_range(1), angle_range(end)]);

    % --- Setup Panel 3: SPR Image ---
    ax3 = nexttile;
    spr_image_matrix = zeros(num_lines, length(angle_range));
    h_spr_image = imagesc(ax3, angle_range, 1:num_lines, spr_image_matrix);
    title(ax3, 'SPR Image', 'FontSize', 14);
    xlabel(ax3, 'Incident Angle (degrees)'); ylabel(ax3, 'Line Index');
    colormap(ax3, 'gray'); caxis(ax3, [0, 1]);
    
    % --- Animation Loop with Variable Speed ---
    slow_step = 20; fast_step = 50;
    interest_window_duration = 250;
    
    fprintf('Generating SPR physics video...\nThis may take some time.\n');
    i = 1; frame_count = 0;
    
    while i <= length(t)
        current_time = t(i);
        is_interesting = any(current_time >= t_breaks & current_time < (t_breaks + interest_window_duration));
        
        if is_interesting, step = slow_step; else, step = fast_step; end
        
        current_rus = line_grams_ru(:, i);
        
        % Update Panel 1 (Sensorgrams)
        for k = 1:num_lines, addpoints(animated_lines(k), t(i), current_rus(k)); end
        time_marker1.Value = t(i);
        
        % Update Panels 2 & 3 (Reflectance and Image)
        for k = 1:num_lines
            % Layer 2 (Analyte) RI is calculated from the RU value
            n2_analyte = n_bulk + (current_rus(k) * RU_TO_RIU);
            
            % n3 is the bulk medium, which is n_bulk
            [Rp_curve, ~] = fresnel_spr_curve(angle_range, n0, n1, n2_analyte, n_bulk, d1, d2, wavelength);
            
            set(reflectance_plots(k), 'YData', Rp_curve);
            spr_image_matrix(k, :) = Rp_curve;
        end
        set(h_spr_image, 'CData', spr_image_matrix);
        
        title(tlo, sprintf('SPR Physics @ t=%.1f s', current_time), 'FontWeight', 'bold', 'FontSize', 18);
        
        drawnow;
        frame = getframe(fig);
        writeVideo(outputVideo, frame);
        
        frame_count = frame_count + 1;
        if mod(frame_count, 10) == 0
            fprintf('  ... generated frame %d (t=%.1f s)\n', frame_count, current_time);
        end
        
        if step < 1; step = 1; end
        i = i + step;
    end
    
    close(outputVideo);
    fprintf('Finished generating %d total frames.\n', frame_count);
    close(fig);
end