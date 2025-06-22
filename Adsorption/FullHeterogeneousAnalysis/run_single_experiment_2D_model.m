function [t_sim, s_sim_matrix] = run_single_experiment_2D_model(p_vector_2D, setting, model_config)
    % This function takes a long parameter vector for all sites, reconstructs
    % the 2D parameter grids, and runs a single simulation.

    % Extract dimensions from model_config
    nx = model_config.gridN_x;
    ny = model_config.gridN_y;
    nz = model_config.gridN_z;
    ads_x = model_config.ads_x_range;
    ads_y = model_config.ads_y_range;
    ads_l = model_config.ads_layer;

    % Calculate the number of sites in the adsorption region
    num_sites = (ads_x(2) - ads_x(1) + 1) * (ads_y(2) - ads_y(1) + 1);

    % --- De-vectorize the parameters ---
    p_kon  = p_vector_2D(1:num_sites);
    p_koff = p_vector_2D(num_sites+1 : 2*num_sites);
    p_smax = p_vector_2D(2*num_sites+1 : 3*num_sites);

    % --- Reconstruct the 2D parameter grids ---
    kon_ads_2D  = reshape(p_kon, [ads_x(2)-ads_x(1)+1, ads_y(2)-ads_y(1)+1]);
    koff_ads_2D = reshape(p_koff, [ads_x(2)-ads_x(1)+1, ads_y(2)-ads_y(1)+1]);
    smax_ads_2D = reshape(p_smax, [ads_x(2)-ads_x(1)+1, ads_y(2)-ads_y(1)+1]);

    % --- Embed the 2D grids into the full 3D simulation grids ---
    kon_grid_full  = zeros(nx, ny, nz);
    koff_grid_full = zeros(nx, ny, nz);
    smax_grid_full = zeros(nx, ny, nz);
    kon_grid_full(ads_x(1):ads_x(2), ads_y(1):ads_y(2), ads_l)  = kon_ads_2D;
    koff_grid_full(ads_x(1):ads_x(2), ads_y(1):ads_y(2), ads_l) = koff_ads_2D;
    smax_grid_full(ads_x(1):ads_x(2), ads_y(1):ads_y(2), ads_l) = smax_ads_2D;

    % --- Run the simulation using the main 3D model ---
    velocity_profile = create_velocity_profile(nz, setting.max_velocity);
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    s0_grid = zeros(nx, ny, nz);

    [t_sim, ~, s_3d_output] = simulate_3d_flow_model_with_pulses(...
        nx, ny, nz, kon_grid_full, koff_grid_full, smax_grid_full, ...
        velocity_profile, t_breaks, concentrations, model_config.D_coeff, ...
        model_config.ru_to_m, s0_grid, model_config.dx, model_config.dz);
    
    % Compute the observable signal (line-averaged)
    s_sim_matrix = compute_s_obs_by_line(s_3d_output, ads_x, ads_y, ads_l);
end