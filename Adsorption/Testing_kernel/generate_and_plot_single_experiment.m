function generate_and_plot_single_experiment()
    % This script generates data for a single in silico experiment,
    % plots the resulting sensorgram with annotations inside the axes, 
    % and saves it as a PNG.
    clearvars; close all; clc;

    % =====================================================================
    % 1. SETUP SIMULATION PARAMETERS
    % =====================================================================
    % Shared parameters
    gridN_x = 22; gridN_y = 5; gridN_z = 3;ads_layer = 1;
    ads_x_range = [5,15]; ads_y_range = [1,5];

    % Fixed physical parameters for Advin
    D_coeff = 6e-3;
    ru_to_m = 1e-6;

    % Experiment timing and concentration settings
    T1 = 2000; T2 = 2*T1; T3 = 3*T1;
    t_total = 4*T1;
    c_diss = 0;
    c1 = 3.3e-6; c2 = 1e-6;

    % Generate settings for ONE specific experiment
    exp_setting = generate_experiments(1, 8.3, T1, T2, T3, t_total, c_diss, c1, c2);
    fprintf('Generated Experiment Settings:\n');
    fprintf('  Concentrations: [%.2e, %.2e, %.2e]\n', exp_setting.pulse_concs);
    fprintf('  Velocity: %.2f\n', exp_setting.max_velocity);


    % =====================================================================
    % 2. GENERATE AND RUN THE EXPERIMENT
    % =====================================================================
    % Create the "ground truth" heterogeneous parameter grids
    [kon_grid_true, koff_grid_true, smax_grid_true] = ...
        create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);

    % Extract only the non-zero adsorption region parameters to pass to the solver
    kon_ads_true = kon_grid_true(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    koff_ads_true = koff_grid_true(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    smax_ads_true = smax_grid_true(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);

    % Run the simulation to get the clean sensorgram data
    [t, s_obs_clean] = run_single_experiment(...
        gridN_x, gridN_y, gridN_z, kon_ads_true, koff_ads_true, smax_ads_true,...
        ads_x_range, ads_y_range, ads_layer, exp_setting, D_coeff, ru_to_m);

    % Add 2% Gaussian noise to create realistic "experimental" data
    noise_level = 0.02;
    rng(101); % for reproducible noise
    s_obs_noisy = s_obs_clean .* (1 + noise_level*randn(size(s_obs_clean)));


    % =====================================================================
    % 3. PLOT AND SAVE THE FIGURE
    % =====================================================================
    figure('Position', [100, 100, 800, 500]);
    hold on;

    % Plot the clean and noisy data
    plot(t, s_obs_noisy, 'b-', 'LineWidth', 2); % Noisy experimental signal
    plot(t, s_obs_clean, 'Color', 'r', 'LineWidth', 0.5); % Faint clean signal

    % --- Add vertical lines and text annotations for each phase ---
    y_lims = ylim;
    % **MODIFIED: Set text position to be 95% of the axis height (inside the plot)**
    text_y_pos = y_lims(1) + 0.95 * (y_lims(2) - y_lims(1)); 

    % Get all time breaks and concentrations for annotation
    all_t_breaks = [0, exp_setting.pulse_times, t_total];
    all_concs = [exp_setting.pulse_concs, c_diss];

    for i = 1:length(all_concs)
        t_start = all_t_breaks(i);
        t_end = all_t_breaks(i+1);

        % Draw vertical line at the end of the phase (except for the last one)
        if i < length(all_concs)
            line([t_end, t_end], y_lims, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);
        end

        % Add text annotation in the middle of the phase
        text_x_pos = (t_start + t_end) / 2;
        text_str = sprintf('C=%.1e M', all_concs(i));
        % **MODIFIED: Added a background to the text for legibility**
        text(text_x_pos, text_y_pos, text_str, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'BackgroundColor', [1 1 1 0.7], ... % White, 70% opaque background
            'Margin', 2); % Padding around text
    end

    % --- Finalize Plot ---
    hold off;
    box on;
    grid on;
    axis tight; % Ensure plot fits data snugly
    %title('Example In Silico Experiment Sensorgram');
    xlabel('Time (s)');
    ylabel('s_{obs}');
    legend('Signal with 2% Noise','Clean Signal','Location', 'best');

    % Save the figure with the correct name
    disp('Saving experiment plot as placeholder_for_experiment_sensorgram.png...');
    print('Testing_kernel/placeholder_for_experiment_sensorgram.png', '-dpng', '-r300');

end


% =========================================================================
% HELPER FUNCTIONS (Copied from the larger script)
% =========================================================================

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

function exp_setting = generate_experiments(M, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2)
    % Simplified to generate just one setting for this script
    exp_setting = struct(...
        'pulse_times', {}, ...
        'pulse_concs', {}, ...
        't_total', {}, ...
        'max_velocity', {}, ...
        'c_diss', {} ...
    );
    exp_setting(1).pulse_times = [T1, T2, T3];
    exp_setting(1).pulse_concs = [c1, 0, c2]; % A standard pattern
    exp_setting(1).t_total = t_total;
    exp_setting(1).max_velocity = base_max_velocity;
    exp_setting(1).c_diss = c_diss;
end

function s_obs = compute_s_obs(s_grid, ads_x_range, ads_y_range, ads_layer)
    ads_cells = s_grid(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    s_obs = squeeze(sum(ads_cells, [2,3,4]));
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

    rng(42); % For reproducibility
    kon_vals = kon_base * (1 + 0.05*randn(ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1));
    koff_vals = koff_base * (1 + 0.05*randn(size(kon_vals)));
    smax_vals = (smax_total/num_ads_cells) * (1 + 0.05*randn(size(kon_vals)));

    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_vals;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_vals;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_vals;
end