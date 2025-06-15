function [t, s_obs_by_line] = run_single_experiment_1D_model(p_1D, setting, model_config)
% Runs a simulation using the 1D heterogeneous parameter model.
    ads_y_dim = model_config.ads_y_range(2) - model_config.ads_y_range(1) + 1;
    ads_x_dim = model_config.ads_x_range(2) - model_config.ads_x_range(1) + 1;
    dx = model_config.dx; dz = model_config.dz;
    % Unpack 1D parameter vector
    kon_1D  = p_1D(1:ads_y_dim);
    koff_1D = p_1D(ads_y_dim+1 : 2*ads_y_dim);
    smax_1D = p_1D(2*ads_y_dim+1 : 3*ads_y_dim);

    % Create 2D parameter grids by replicating parameters along the x-axis
    kon_ads  = repmat(kon_1D', ads_x_dim, 1);
    koff_ads = repmat(koff_1D', ads_x_dim, 1);
    smax_ads = repmat(smax_1D', ads_x_dim, 1);
    
    % Create full 3D parameter grids
    [kon_grid, koff_grid, smax_grid] = create_heterogeneous_grids_from_ads(...
        model_config.gridN_x, model_config.gridN_y, model_config.gridN_z, ...
        model_config.ads_x_range, model_config.ads_y_range, model_config.ads_layer, ...
        kon_ads, koff_ads, smax_ads);

    % Create velocity, time, and concentration profiles
    velocity_profile = create_velocity_profile(model_config.gridN_z, setting.max_velocity);
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    s0_grid = zeros(model_config.gridN_x, model_config.gridN_y, model_config.gridN_z);

    % Run simulation
    [t, ~, s] = simulate_3d_flow_model_with_pulses(...
        model_config.gridN_x, model_config.gridN_y, model_config.gridN_z, ...
        kon_grid, koff_grid, smax_grid, velocity_profile, ...
        t_breaks, concentrations, model_config.D_coeff, model_config.ru_to_m, s0_grid,dx,dz);
    
    % Compute the per-line observed signal
    s_obs_by_line = compute_s_obs_by_line(s, model_config.ads_x_range, ...
        model_config.ads_y_range, model_config.ads_layer);
end