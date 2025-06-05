function validate_kernel_implementation_with_experiments()
    clearvars; close all; clc;
    % Shared parameters
    gridN_x = 25; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [10,14]; ads_y_range = [2,4];
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Homogeneous parameters (from step 1 results)
    homog_params = [9.4e3, 0.0078, 1.0]; % kon, koff, smax_total
    smax_per_cell = homog_params(3) / num_ads_cells;
    
    % Fixed parameters
    D_coeff = 6e-3;
    ru_to_m = 1e-6;
    base_max_velocity = 8.3;
    
    % Example pulse parameters
    T1 = 800; T2 = 1200; T3 = 2000;
    c1 = 3.3e-6; c2 = 1e-6; c_diss = 0;
    
    % Experiment settings
    exp_settings = struct(...
        'pulse_times', { [T1, T2, T3], [T1, T2, T3], [T1, T2, T3] }, ...
        'pulse_concs', { [c1, 0, c2], [c1, 0, 2*c2], [2*c1, 0, c2] }, ...
        't_total', {2400, 2400, 2400}, ...
        'max_velocity', {base_max_velocity, 2*base_max_velocity, 0.5*base_max_velocity}, ...
        'c_diss', {c_diss, c_diss, c_diss} ...
    );
    
    % Create heterogeneous parameters (same for all experiments)
    [kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
        create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
    for exp_idx = 1:3
        setting = exp_settings(exp_idx);
        fprintf('\nRunning Experiment %d\n', exp_idx);
        
        % Create homogeneous grids
        [kon_grid_homog, koff_grid_homog, smax_grid_homog] = ...
            create_homogeneous_grids(homog_params, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
        
        % Create velocity profile
        velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);
        
        % Time breaks and concentrations
        t_breaks = [0, setting.pulse_times, setting.t_total];
        concentrations = [setting.pulse_concs, setting.c_diss];
        
        % Initial condition for s (empty)
        s0_grid = zeros(gridN_x, gridN_y, gridN_z);
        
        % Homogeneous simulation
        [t_homog, ~, s_homog, ~, ~] = simulate_3d_flow_model_with_pulses(...
            gridN_x, gridN_y, gridN_z, kon_grid_homog, koff_grid_homog, smax_grid_homog, ...
            velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid);
        
        % Heterogeneous simulation
        [t_heterog, ~, s_heterog, ~, ~] = simulate_3d_flow_model_with_pulses(...
            gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
            velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid);
        
        % Compute s_obs (sum over adsorption region)
        s_obs_homog = compute_s_obs(s_homog, ads_x_range, ads_y_range, ads_layer);
        s_obs_heterog = compute_s_obs(s_heterog, ads_x_range, ads_y_range, ads_layer);
        
        % Plot results
        plot_composite_behavior(t_homog, s_obs_homog, t_heterog, s_obs_heterog, ...
            setting.pulse_times, setting.max_velocity, exp_idx);
    end
end

%% New Helper Functions
function s_obs = compute_s_obs(s_grid, ads_x_range, ads_y_range, ads_layer)
    ads_cells = s_grid(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    s_obs = squeeze(sum(ads_cells, [2,3,4]));
end

function plot_composite_behavior(t_homog, s_obs_homog, t_heterog, s_obs_heterog, pulse_times, max_velocity, exp_idx)
    figure('Position', [100, 100, 800, 600]);
    plot(t_homog, s_obs_homog, 'b-', 'LineWidth', 2); hold on;
    plot(t_heterog, s_obs_heterog, 'r--', 'LineWidth', 1.5);
    
    % Mark pulse transitions
    y_lims = [min([s_obs_homog; s_obs_heterog]), max([s_obs_homog; s_obs_heterog])];
    for i = 1:length(pulse_times)
        line([pulse_times(i), pulse_times(i)], y_lims, ...
            'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
    end
    
    xlabel('Time (s)'); 
    ylabel('s_{obs}(t)');
    title(sprintf('Exp %d: Composite Behavior (v_{max}=%.1f)', exp_idx, max_velocity));
    legend('Homogeneous', 'Heterogeneous (5% var)', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);
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
%             'pulse_times', { [T1, T2, T3], [T1, T2, T3], [T1, T2, T3] }, ...        t_breaks = [0, setting.pulse_times, setting.t_total];


    % Process each time segment
    num_segments = length(t_breaks) - 1;
    for seg = 1:num_segments
        t_start = t_breaks(seg);
        t_end = t_breaks(seg+1);
        c0_seg = concentrations(seg);
        disp(t_start);disp(c0_seg)
        % Determine time points for segment (min 10 points)
        num_points = max(10, ceil(200 * (t_end - t_start) / (t_breaks(end) - t_breaks(1))));
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
        if seg<num_segments
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
    % Reshape state variables
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);       % Corrected indices for s
    Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);    % Corrected indices for Q
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
    
    % Adsorption kinetics with surface capacity constraint
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
%% Helper Functions
function [kon_grid, koff_grid, smax_grid] = create_homogeneous_grids(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    % Create uniform parameter grids
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
    
    % Initialize grids
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    % Add 5% variation to adsorption region
    rng(42); % For reproducibility
    kon_vals = kon_base * (1 + 0.05*randn(ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1));
    koff_vals = koff_base * (1 + 0.05*randn(size(kon_vals)));
    smax_vals = (smax_total/num_ads_cells) * (1 + 0.05*randn(size(kon_vals)));
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_vals;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_vals;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_vals;
end