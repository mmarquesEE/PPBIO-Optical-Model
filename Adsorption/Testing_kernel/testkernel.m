%% KERNEL VALIDATION SCRIPT
function validate_kernel_implementation()
    clearvars; close all; clc;
    % Shared parameters
    gridN_x = 25; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [10,14]; ads_y_range = [2,4];
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Homogeneous parameters
    homog_params = [9.4e3, 0.0078, 1.0]; % kon, koff, smax_total
    smax_per_cell = homog_params(3) / num_ads_cells;
    
    % 1. Homogeneous Kernel Calculation
    [K_homog, t_homog, s_homog, Q_homog] = simulate_homogeneous_kernel(homog_params, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
    % 2. Heterogeneous Kernel Calculation
    [kon_grid_het, koff_grid_het, smax_grid_het] = create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    [K_heterog, t_heterog, s_heterog, Q_heterog] = simulate_heterogeneous_kernel(kon_grid_het, koff_grid_het, smax_grid_het, gridN_x, gridN_y, gridN_z);
    
    % --- Generate all grids for plotting ---
    [kon_grid_hom, koff_grid_hom, smax_grid_hom] = create_homogeneous_grids(homog_params, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
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
end

%% Helper Functions

% =========================================================================
% UPDATED FUNCTION TO PLOT AND SAVE ALL PARAMETER FIELDS
% =========================================================================
function plot_and_save_parameter_fields(kon_homog, koff_homog, smax_homog, ...
                                        kon_heterog, koff_heterog, smax_heterog, ...
                                        ads_x_range, ads_y_range, ads_layer)
    
    figure('Position', [300, 300, 500, 400]);
    
    % Use tiledlayout for better control over spacing
    tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Row 1: Plot k_on ---
    % Calculate color limits ONLY from the sensible region
    kon_slice_heterog = kon_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_kon = [min(kon_slice_heterog(:)), max(kon_slice_heterog(:))];
    if diff(c_limits_kon) < 1e-9; c_limits_kon(2) = c_limits_kon(1) + 1; end % Robustness check

    % Homogeneous k_on
    nexttile;
    plot_grid_with_black_background(kon_homog(:,:,ads_layer)', c_limits_kon);
    title('Homogeneous');
    ylabel('k_{on}', 'FontSize', 10);
    set(gca, 'XTickLabel', []); 
    
    % Heterogeneous k_on
    nexttile;
    plot_grid_with_black_background(kon_heterog(:,:,ads_layer)', c_limits_kon);
    title('Heterogeneous');
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    colorbar;

    % --- Row 2: Plot k_off ---
    koff_slice_heterog = koff_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_koff = [min(koff_slice_heterog(:)), max(koff_slice_heterog(:))];
    if diff(c_limits_koff) < 1e-9; c_limits_koff(2) = c_limits_koff(1) + 1; end

    % Homogeneous k_off
    nexttile;
    plot_grid_with_black_background(koff_homog(:,:,ads_layer)', c_limits_koff);
    ylabel('k_{off}', 'FontSize', 10);
    set(gca, 'XTickLabel', []); 

    % Heterogeneous k_off
    nexttile;
    plot_grid_with_black_background(koff_heterog(:,:,ads_layer)', c_limits_koff);
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    colorbar;

    % --- Row 3: Plot s_max ---
    smax_slice_heterog = smax_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_smax = [min(smax_slice_heterog(:)), max(smax_slice_heterog(:))];
    if diff(c_limits_smax) < 1e-9; c_limits_smax(2) = c_limits_smax(1) + 1; end

    % Homogeneous s_max
    nexttile;
    plot_grid_with_black_background(smax_homog(:,:,ads_layer)', c_limits_smax);
    xlabel('y-grid index');
    ylabel('s_{max}', 'FontSize', 10);
    
    % Heterogeneous s_max
    nexttile;
    plot_grid_with_black_background(smax_heterog(:,:,ads_layer)', c_limits_smax);
    xlabel('y-grid index');
    set(gca, 'YTickLabel', []);
    colorbar;
    
    % Add a main y-label for the entire layout
    han = gcf();
    han.CurrentAxes = gca();
    ylabel(han.CurrentAxes.Parent, 'x-grid index', 'FontSize',10)
    
    % Save the figure
    disp('Saving parameter fields plot as param_fields.png...');
    print('Testing_kernel/param_fields.png', '-dpng', '-r300');
end

% --- Helper sub-function for plotting ---
function plot_grid_with_black_background(grid_data, c_limits)
    % This function uses transparency to make zero-value areas reveal a black background
    
    % Plot the image and get a handle to it
    h = imagesc(grid_data);
    
    % Set the colormap for the data
    colormap(gca, parula);
    
    % Create a transparency map: 1 for non-zero data, 0 for zero-data
    alpha_map = double(grid_data ~= 0);
    
    % Apply the transparency map
    set(h, 'AlphaData', alpha_map);
    
    % Set the axis background color to black
    set(gca, 'Color', 'k');
    
    % Apply the color limits and tighten the axis
    caxis(c_limits);
    axis tight;
end


function [kon_grid, koff_grid, smax_grid] = create_homogeneous_grids(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
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

function [K, t, s, Q] = simulate_homogeneous_kernel(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    [kon_grid, koff_grid, smax_grid] = create_homogeneous_grids(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer);
    velocity_profile = create_velocity_profile(nz, 8.3);
    
    % Create initial condition for s (20% saturation in adsorption region)
    s0_grid = 0.2 * smax_grid;
    
    [t, ~, s, K, Q] = simulate_3d_flow_model(...
        nx, ny, nz, kon_grid, koff_grid, smax_grid,...
        velocity_profile, 3.3e-6, 0, 1500, 3500, 6e-3, 1e-6, s0_grid);
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

function [K, t, s, Q] = simulate_heterogeneous_kernel(kon_grid, koff_grid, smax_grid, nx, ny, nz)
    velocity_profile = create_velocity_profile(nz, 8.3);
    
    % Create initial condition for s (20% saturation in adsorption region)
    s0_grid = 0.2 * smax_grid;
    
    [t, ~, s, K, Q] = simulate_3d_flow_model(...
        nx, ny, nz, kon_grid, koff_grid, smax_grid,...
        velocity_profile, 3.3e-6, 0, 1500, 3500, 6e-3, 1e-6, s0_grid);
end

% =========================================================================
% MODIFIED COMPARISON FUNCTION TO SAVE PLOTS
% =========================================================================
function compare_composite_behavior(s_homog, s_heterog, Q_homog, Q_heterog, K_homog, K_heterog,...
        t_homog, t_heterog, ads_x_range, ads_y_range, ads_layer,...
        kon_homog, smax_per_cell, kon_heterog_ads, smax_heterog_ads)
    
    % Extract adsorption region for homogeneous
    ads_cells_s_homog = s_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_Q_homog = Q_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_K_homog = K_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    % Extract adsorption region for heterogeneous
    ads_cells_s_heterog = s_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_Q_heterog = Q_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_K_heterog = K_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    s0_homog_ads = squeeze(ads_cells_s_homog(1, :, :));
    s0_heterog_ads = squeeze(ads_cells_s_heterog(1, :, :));
    
    alpha_homog_ads = kon_homog * smax_per_cell * ones(size(kon_heterog_ads));
    alpha_heterog_ads = kon_heterog_ads .* smax_heterog_ads;
    
    alpha_homog_ads = reshape(alpha_homog_ads, [1, size(alpha_homog_ads)]);
    alpha_heterog_ads = reshape(alpha_heterog_ads, [1, size(alpha_heterog_ads)]);
    s0_homog_ads = reshape(s0_homog_ads, [1, size(s0_homog_ads)]);
    s0_heterog_ads = reshape(s0_heterog_ads, [1, size(s0_heterog_ads)]);

    decay_term_homog = s0_homog_ads .* exp(-ads_cells_Q_homog);
    decay_term_heterog = s0_heterog_ads .* exp(-ads_cells_Q_heterog);
    
    kernel_term_homog = alpha_homog_ads .* ads_cells_K_homog;
    kernel_term_heterog = alpha_heterog_ads .* ads_cells_K_heterog;
    
    decay_sum_homog = squeeze(sum(decay_term_homog, [2,3]));
    kernel_sum_homog = squeeze(sum(kernel_term_homog, [2,3]));
    s_obs_homog = decay_sum_homog + kernel_sum_homog;
    
    decay_sum_heterog = squeeze(sum(decay_term_heterog, [2,3]));
    kernel_sum_heterog = squeeze(sum(kernel_term_heterog, [2,3]));
    s_obs_heterog = decay_sum_heterog + kernel_sum_heterog;
    
    discrepancy = trapz(t_heterog, (s_obs_homog - s_obs_heterog).^2);
    fprintf('Composite behavior discrepancy: %.2e\n', discrepancy);
    
    % Plot decomposition for both cases
    fig1 = figure;
    subplot(2,1,1);
    plot(t_homog, s_obs_homog, 'k-', 'LineWidth', 2); hold on;
    plot(t_homog, decay_sum_homog, 'b--', 'LineWidth', 1.5);
    plot(t_homog, kernel_sum_homog, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('s_{obs}');
    legend('Total', 'Decay Term', 'Kernel Term', 'Location', 'best');
    %title('(A) Homogeneous Case: Signal Decomposition');
    title('(A)')
    grid on;
    
    subplot(2,1,2);
    plot(t_heterog, s_obs_heterog, 'k-', 'LineWidth', 2); hold on;
    plot(t_heterog, decay_sum_heterog, 'b--', 'LineWidth', 1.5);
    plot(t_heterog, kernel_sum_heterog, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('s_{obs}');
    legend('Total', 'Decay Term', 'Kernel Term', 'Location', 'best');
    %title('(B) Heterogeneous Case: Signal Decomposition');
    title('(B)')
    grid on;
    
    disp('Saving signal decomposition plot as signal_decomposition.png...');
    print(fig1, 'Testing_kernel/signal_decomposition.png', '-dpng', '-r300');
    
    % Plot total signal comparison
    fig2 = figure;
    plot(t_heterog, s_obs_homog, 'b-', 'LineWidth', 2); hold on;
    plot(t_heterog, s_obs_heterog, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('s_{obs}(t)');
    legend('Homogeneous', 'Heterogeneous (5% var)', 'Location', 'best');
    %title(sprintf('Composite Behavior Comparison\nDiscrepancy: %.2e', discrepancy));
    grid on;
    
    disp('Saving sensorgram comparison plot as sensorgram_comparison.png...');
    print(fig2, 'Testing_kernel/sensorgram_comparison.png', '-dpng', '-r300');
end

function [t, c_s, s, K, Q] = simulate_3d_flow_model(nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff, ru_to_m, s0_grid)
    % Initialize concentrations and auxiliary variables
    c_s = zeros(nx, ny, nz);
    c_s(1, :, :) = c0_assoc;
    s = s0_grid;  % Use provided initial condition
    Q = zeros(nx, ny, nz);
    R = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:); Q(:); R(:)];

    % Time parameters
    tspan_assoc = linspace(0, t_assoc, 800);
    tspan_diss = linspace(t_assoc, t_total, 800);
    
    % Solve ODE for association phase
    options = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0_assoc, D_coeff, ru_to_m), tspan_assoc, y0, options);
    
    % Reset for dissociation phase
    y_end_assoc = y_assoc(end,:)';
    num_cells = nx*ny*nz;
    c_s_end = reshape(y_end_assoc(1:num_cells), [nx, ny, nz]);
    s_end = reshape(y_end_assoc(num_cells+1:2*num_cells), [nx, ny, nz]);
    c_s_end(1, :, :) = c0_diss;
    y0_diss = [c_s_end(:); s_end(:); y_end_assoc(2*num_cells+1:end)]; % Preserve Q/R
    
    % Solve ODE for dissociation phase
    [t_diss, y_diss] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0_diss, D_coeff, ru_to_m), tspan_diss, y0_diss, options);
    
    % Combine results
    t = [t_assoc; t_diss(2:end)];
    y = [y_assoc; y_diss(2:end,:)];
    
    % Extract variables
    c_s = reshape(y(:,1:num_cells), [length(t), nx, ny, nz]);
    s = reshape(y(:,num_cells+1:2*num_cells), [length(t), nx, ny, nz]);
    Q = reshape(y(:,2*num_cells+1:3*num_cells), [length(t), nx, ny, nz]);
    R = reshape(y(:,3*num_cells+1:end), [length(t), nx, ny, nz]);
    
    % Compute kernel K(x,y,t) = exp(-Q) .* R
    K = exp(-Q) .* R;
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