clearvars; close all; clc;
    % Shared parameters
gridN_x = 22; gridN_y = 5; gridN_z = 3;ads_layer = 1;
ads_x_range = [5,15]; ads_y_range = [1,5];
num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);

% Homogeneous parameters
homog_params = [9.4e3, 0.0078, 2960]; % kon, koff, smax_total
smax_per_cell = homog_params(3) / num_ads_cells;
velocity_profile = create_velocity_profile(gridN_z, 8.3);
% Fixed parameters
D_coeff = 6e-5;ru_to_m = 1e-10;
grid_size_x = 11.0; % mm
grid_size_z = 0.3; % mm
dx = grid_size_x / gridN_x; % Size of one grid cell in x-direction (mm)
dz = grid_size_z / gridN_z; % Size of one grid cell in z-direction (mm)
T1 = 1000; T2 = 2*T1;
c_diss = 0;c1 = 3.3e-6;

% 1. Homogeneous Kernel Calculation
[kon_grid_hom, koff_grid_hom, smax_grid_hom] = create_homogeneous_grids(homog_params, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
s0_grid_hom = 0.2 * smax_grid_hom;
[t_homog, ~, s_homog, K_homog, Q_homog] = simulate_3d_flow_model_ALL(...
    gridN_x, gridN_y, gridN_z, kon_grid_hom, koff_grid_hom, smax_grid_hom,...
    velocity_profile,c1, c_diss, T1, T2, D_coeff, ru_to_m, s0_grid_hom,dx,dz);
% 2. Heterogeneous Kernel Calculation
[kon_grid_het, koff_grid_het, smax_grid_het] = create_ground_truth_heterogeneity_full(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
s0_grid_het = 0.2 * smax_grid_het;
[t_heterog, ~, s_heterog, K_heterog, Q_heterog] = simulate_3d_flow_model_ALL(...
gridN_x, gridN_y, gridN_z, kon_grid_het, koff_grid_het, smax_grid_het,...
velocity_profile,c1, c_diss, T1, T2, D_coeff, ru_to_m, s0_grid_het,dx,dz);    
% --- Plot all Parameter Fields (Updated) ---
plot_and_save_parameter_fields(kon_grid_hom, koff_grid_hom, smax_grid_hom, ...
                               kon_grid_het, koff_grid_het, smax_grid_het, ...
                               ads_x_range, ads_y_range, ads_layer);

% Extract adsorption region parameters for comparison function
kon_heterog_ads = kon_grid_het(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
smax_heterog_ads = smax_grid_het(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);

% 3. Composite Behavior Comparison and Plotting
compare_composite_behavior(s_homog, s_heterog, Q_homog, Q_heterog, K_homog, K_heterog,...
    t_homog, t_heterog, ads_x_range, ads_y_range, ads_layer,...
    homog_params(1), smax_per_cell, kon_heterog_ads, smax_heterog_ads);