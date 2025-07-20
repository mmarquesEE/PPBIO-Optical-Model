%% Main Script to Run and Generate a Professional-Grade, Multi-Panel Video
% This script runs the 3D simulation and generates a final, high-resolution,
% high-framerate MP4 video with a custom light gray background.
%
clearvars; close all; clc;

fprintf('--- Setting up Professional-Grade Hierarchical Video Animation ---\n');

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

% Get the data in [nx, ny, nz, time] format
[t, C_grid_t, s_grid_t] = simulate_3d_flow_model_with_pulses_fluid(...
    gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
    velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid, dx, dz);

fprintf('Simulation complete. Generated %d time points.\n', length(t));
toc;

%% --- STEP 3: Pre-process Data for Hierarchical Animation ---
% ========================================================================
fprintf('Pre-processing data for animation...\n');

s_ads_region_t = squeeze(s_grid_t(ads_x_range(1):ads_x_range(2), :, ads_layer, :));
line_grams = squeeze(sum(s_ads_region_t, 1));
line_to_inspect = round(gridN_y / 2);
site_grams = squeeze(s_ads_region_t(:, line_to_inspect, :));
C_grid_for_anim = permute(C_grid_t, [4, 1, 2, 3]);
s_surface_for_anim = squeeze(s_grid_t(:,:,ads_layer,:));

%% --- STEP 4: Call the Comprehensive Animation Function ---
% ========================================================================
create_hierarchical_animation(t, C_grid_for_anim, s_surface_for_anim, ...
                              line_grams, site_grams, ...
                              gridN_x, gridN_y, gridN_z, t_breaks, concentrations, line_to_inspect);

fprintf('\nProfessional-grade video animation saved successfully.\n');


%% ========================================================================
%  --- LOCAL FUNCTION: Hierarchical Multi-Panel Animation (Pro Video) ---
%  ========================================================================
function create_hierarchical_animation(t, c_s, s_surf, line_grams, site_grams, gridN_x, gridN_y, gridN_z, t_breaks, concentrations, line_to_inspect)
    
    % --- Setup VideoWriter for high-quality MP4 output ---
    filename = 'hierarchical_simulation_pro.mp4';
    outputVideo = VideoWriter(filename, 'MPEG-4');
    outputVideo.Quality = 100; % Maximum quality
    outputVideo.FrameRate = 30; % Cinema-smooth frame rate
    open(outputVideo);

    % --- NEW: Set the custom light gray background color ---
    background_color = [217/255, 217/255, 217/255]; % Light gray (#D9D9D9)
    
    % --- Increased figure size for higher resolution ---
    fig = figure('Position', [50 50 1920 1440], 'Color', background_color);
    
    tlo = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Setup Panel 1: 3D Fluid Flow ---
    ax1 = nexttile;
    max_conc = max(c_s(:)); if max_conc == 0; max_conc = 1; end
    [X,Y,Z] = meshgrid(1:gridN_x, 1:gridN_y, 1:gridN_z);
    sx = [1, gridN_x]; sy = round(gridN_y/2); sz = [1, gridN_z];
    conc3d_permuted = permute(squeeze(c_s(1,:,:,:)), [2 1 3]);
    h_slice = slice(ax1, X,Y,Z, conc3d_permuted, sx, sy, sz);
    set(h_slice, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    hold(ax1, 'on');
    [inX, inY, inZ] = meshgrid(1, 1:gridN_y, 1:gridN_z);
    scatter3(ax1, inX(:), inY(:), inZ(:), 30, 'g', 'filled', 'MarkerFaceAlpha', 0.6);
    [outX, outY, outZ] = meshgrid(gridN_x, 1:gridN_y, 1:gridN_z);
    scatter3(ax1, outX(:), outY(:), outZ(:), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    hold(ax1, 'off');
    shading(ax1, 'interp'); colormap(ax1, 'parula'); caxis(ax1, [0 max_conc]);
    title(ax1, '3D Analyte Concentration', 'FontSize', 14);
    xlabel(ax1, 'Length (x)'); ylabel(ax1, 'Depth (y)'); zlabel(ax1, 'Height (z)');
    view(ax1, 45, 25); axis(ax1, 'tight');
    set(ax1, 'Color', 'none'); % Make axes background transparent to the figure color

    % --- Setup Panel 2: 2D Surface Coverage ---
    ax2 = nexttile;
    h_surf = imagesc(ax2, s_surf(:,:,1)');
    title(ax2, 'Sensor Surface Coverage (Top View)', 'FontSize', 14);
    xlabel(ax2, 'Length (x)'); ylabel(ax2, 'Line Index (y)');
    colormap(ax2, 'hot'); caxis(ax2, [0, max(s_surf(:)) + eps]); colorbar(ax2);
    axis(ax2, 'xy');
    set(ax2, 'Color', 'none');

    % --- Setup Panel 3: Line Sensorgrams ---
    ax3 = nexttile;
    hold(ax3, 'on');
    line_colors = lines(gridN_y);
    animated_lines = gobjects(gridN_y, 1);
    for k = 1:gridN_y, animated_lines(k) = animatedline(ax3, 'Color', line_colors(k,:), 'LineWidth', 2); end
    hold(ax3, 'off'); grid(ax3, 'on');
    title(ax3, 'Line-Level Sensorgrams', 'FontSize', 14); xlabel(ax3, 'Time (s)'); ylabel(ax3, 'Line Response');
    xlim(ax3, [0, t(end)]);
    legend(ax3, arrayfun(@(k) sprintf('Line %d', k), 1:gridN_y, 'UniformOutput', false), 'Location', 'northwest');
    time_marker3 = xline(ax3, 0, 'r--', 'LineWidth', 1.5);
    set(ax3, 'Color', 'none');

    % --- Setup Panel 4: Site Sensorgrams for a single line ---
    ax4 = nexttile;
    num_sites = size(site_grams, 1);
    site_colors = jet(num_sites);
    animated_sites = gobjects(num_sites, 1);
    hold(ax4, 'on');
    for k = 1:num_sites, animated_sites(k) = animatedline(ax4, 'Color', site_colors(k,:), 'LineWidth', 1.5); end
    hold(ax4, 'off'); grid(ax4, 'on');
    title(ax4, sprintf('Site-Level Sensorgrams for Line %d', line_to_inspect), 'FontSize', 14);
    xlabel(ax4, 'Time (s)'); ylabel(ax4, 'Site Response');
    xlim(ax4, [0, t(end)]);
    time_marker4 = xline(ax4, 0, 'r--', 'LineWidth', 1.5);
    set(ax4, 'Color', 'none');

    % --- Animation Loop with Variable Speed & Annotations ---
    slow_step = 20; fast_step = 50;
    interest_window_duration = 250; % Extended slow-motion duration
    
    fprintf('Generating professional-quality video...\nThis may take some time.\n');
    i = 1; frame_count = 0;
    phase_text_handle = [];
    
    while i <= length(t)
        current_time = t(i);
        is_interesting = any(current_time >= t_breaks & current_time < (t_breaks + interest_window_duration));
        
        if is_interesting, step = slow_step; else, step = fast_step; end
        
        current_phase_idx = find(t_breaks <= current_time, 1, 'last');
        if i > step, previous_phase_idx = find(t_breaks <= t(i-step), 1, 'last');
        else, previous_phase_idx = 0; end

        if current_phase_idx > previous_phase_idx
            delete(phase_text_handle);
            current_conc = concentrations(current_phase_idx);
            if current_conc > 0, phase_str = sprintf('Injection: C = %.1e M', current_conc); text_color = 'blue';
            else, phase_str = 'Washout Phase'; text_color = 'black'; end
            phase_text_handle = text(ax1, 0.05, 0.9, phase_str, 'Units', 'normalized', ...
                'FontSize', 16, 'FontWeight', 'bold', 'Color', text_color, 'BackgroundColor', [1 1 1 0.7]);
        end
        
        % Update all panels
        delete(h_slice);
        conc3d_permuted = permute(squeeze(c_s(i,:,:,:)), [2 1 3]);
        h_slice = slice(ax1, X,Y,Z, conc3d_permuted, sx, sy, sz);
        set(h_slice, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        set(h_surf, 'CData', s_surf(:,:,i)');
        for k = 1:gridN_y, addpoints(animated_lines(k), t(i), line_grams(k, i)); end
        for k = 1:num_sites, addpoints(animated_sites(k), t(i), site_grams(k, i)); end
        time_marker3.Value = t(i);
        time_marker4.Value = t(i);
        title(tlo, sprintf('Simulation Dashboard @ t=%.1f s', current_time), 'FontWeight', 'bold', 'FontSize', 18);
        
        % Write frame to video file
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
    
    % Finalize and close video file
    close(outputVideo);
    fprintf('Finished generating %d total frames.\n', frame_count);
    close(fig);
end