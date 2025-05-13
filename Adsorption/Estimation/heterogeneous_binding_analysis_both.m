function heterogeneous_binding_analysis()
    % Main function for analyzing heterogeneous binding parameters
    clearvars; close all; clc;
    
    % ================ SIMULATION PARAMETERS ================
    gridN_x = 25;          % Reduced grid for faster testing
    gridN_y = 5;
    gridN_z = 3;
    ads_layer = 1;         % Adsorption layer (z=1)
    ads_x_range = [10,15]; % Central adsorption region (6 cells)
    ads_y_range = [2,4];   % 3 cells in y-direction
    c0_assoc = 3.3e-6;     % Association phase concentration
    c0_diss = 0;           % Dissociation phase concentration
    t_association = 1500;  % Reduced time for testing
    t_dissociation = 2000;
    D_coeff = 6e-3;        % Diffusion coefficient
    ru_to_m = 1e-6;        % RU conversion factor
    
    % Flow parameters
    velocity_profile = create_velocity_profile(gridN_z, 8.3);
    
    % ================ HETEROGENEOUS PARAMETER SETUP ================
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells) * length(ads_y_cells);
    
    % ================ SYNTHETIC DATA GENERATION ================
    fprintf('Generating synthetic data...\n');
    rng(1);  % For reproducibility
    
    % True parameters (log-uniform distributions)
    true_kon = 10.^(rand(num_ads_cells,1)*2 - 1);  % 1e-1 to 1e4
    true_koff = 10.^(rand(num_ads_cells,1)*2 - 2); % 1e-5 to 1e-2
    true_kon_grid = reshape(true_kon, length(ads_x_cells), length(ads_y_cells));
    true_koff_grid = reshape(true_koff, length(ads_x_cells), length(ads_y_cells));
    
    % Create parameter grids
    [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        true_kon_grid, true_koff_grid);
    
    % Run forward simulation
    [t, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_grid,...
        koff_grid, smax_grid, velocity_profile, c0_assoc, c0_diss,...
        t_association, t_dissociation, D_coeff, ru_to_m);
    
    % Add noise to observed signal
    s_obs = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    noise_level = 0.02;
    s_obs = s_obs + noise_level*max(s_obs)*randn(size(s_obs));
    
    % ================ INVERSION SETUP ================
    fprintf('Setting up inverse problem...\n');
    
    % Combined parameter vector (log-scale)
    p0_kon = log10(true_kon_grid(:)) + 0.1*randn(num_ads_cells,1);
    p0_koff = log10(true_koff_grid(:)) + 0.1*randn(num_ads_cells,1);
    p0 = [p0_kon; p0_koff];
    
    % Parameter bounds (log-scale)
    lb_kon = log10(1e-2 * ones(num_ads_cells,1));
    ub_kon = log10(1e5 * ones(num_ads_cells,1));
    lb_koff = log10(1e-5 * ones(num_ads_cells,1));
    ub_koff = log10(1e0 * ones(num_ads_cells,1));
    lb = [lb_kon; lb_koff];
    ub = [ub_kon; ub_koff];
    
    % Regularization parameters
    lambda_kon = 1e-3;
    lambda_koff = 1e-4;

    % ================ SOLVE INVERSION ================
    fprintf('Solving inverse problem...\n');
    options = optimoptions('fmincon', 'Display', 'iter',...
        'Algorithm', 'interior-point', 'MaxIterations', 30,...
        'UseParallel', true);
    
    p_opt = fmincon(@(p) cost_function(p, s_obs, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, c0_assoc,...
        c0_diss, t_association, t_dissociation, D_coeff, ru_to_m,...
        lambda_kon, lambda_koff), p0, [], [], [], [], lb, ub, [], options);
    
    % Split optimized parameters
    recovered_kon = 10.^p_opt(1:num_ads_cells);
    recovered_koff = 10.^p_opt(num_ads_cells+1:end);
    recovered_kon_grid = reshape(recovered_kon, size(true_kon_grid));
    recovered_koff_grid = reshape(recovered_koff, size(true_koff_grid));

    % ================ VISUALIZATION ================
    % Parameter maps
    figure;
    subplot(2,2,1);
    imagesc(log10(true_kon_grid'));
    colorbar; title('True log(k_{on})'); axis equal tight;
    
    subplot(2,2,2);
    imagesc(log10(recovered_kon_grid'));
    colorbar; title('Recovered log(k_{on})'); axis equal tight;
    
    subplot(2,2,3);
    imagesc(log10(true_koff_grid'));
    colorbar; title('True log(k_{off})'); axis equal tight;
    
    subplot(2,2,4);
    imagesc(log10(recovered_koff_grid'));
    colorbar; title('Recovered log(k_{off})'); axis equal tight;

    % Sensorgram comparison
    [~, s_fit] = cost_function(p_opt, s_obs, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, c0_assoc,...
        c0_diss, t_association, t_dissociation, D_coeff, ru_to_m, 0, 0);
    
    figure;
    plot(t, s_obs, 'o', 'MarkerSize', 4, 'DisplayName', 'Observed');
    hold on;
    plot(t, s_fit, 'LineWidth', 2, 'DisplayName', 'Fit');
    xlabel('Time (s)'); ylabel('Response (RU)');
    legend; title('Sensorgram Comparison');
end

%% Helper Functions
function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz - 1);
    h = nz - 1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1, 1, nz]);
end

function [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_values, koff_values)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    x_idx = ads_x_range(1):ads_x_range(2);
    y_idx = ads_y_range(1):ads_y_range(2);
    
    kon_grid(x_idx, y_idx, ads_layer) = kon_values;
    koff_grid(x_idx, y_idx, ads_layer) = koff_values;
    smax_grid(x_idx, y_idx, ads_layer) = 1;  % Fixed surface capacity
end

function [Dx, Dy] = create_tv_operators(nx, ny)
    % Horizontal differences
    Dx = spdiags([-ones(nx,1) ones(nx,1)], [0 1], nx, nx);
    Dx = kron(speye(ny), Dx);
    
    % Vertical differences
    Dy = spdiags([-ones(ny,1) ones(ny,1)], [0 1], ny, ny);
    Dy = kron(Dy, speye(nx));
end

%% Cost Function with Regularization
function [cost, s_fit] = cost_function(p, s_obs, nx, ny, nz, ads_x_range,...
    ads_y_range, ads_layer, velocity_profile, c0_assoc, c0_diss,...
    t_assoc, t_total, D_coeff, ru_to_m, lambda_kon, lambda_koff)
    
    % Split parameters
    num_ads_cells = length(ads_x_range(1):ads_x_range(2)) * length(ads_y_range(1):ads_y_range(2));
    kon_params = p(1:num_ads_cells);
    koff_params = p(num_ads_cells+1:end);
    
    % Reshape to grids
    num_x = length(ads_x_range(1):ads_x_range(2));
    num_y = length(ads_y_range(1):ads_y_range(2));
    kon_grid = reshape(10.^kon_params, num_x, num_y);
    koff_grid = reshape(10.^koff_params, num_x, num_y);
    
    % Generate 3D parameter grids
    [kon_3d, koff_3d, smax_3d] = generate_hetero_param_grids(...
        nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_grid, koff_grid);
    
    % Run forward model
    [~, ~, s] = simulate_3d_flow_model(nx, ny, nz, kon_3d, koff_3d,...
        smax_3d, velocity_profile, c0_assoc, c0_diss, t_assoc,...
        t_total, D_coeff, ru_to_m);
    
    % Extract signal
    s_fit = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    
    % Regularization terms
    [Dx, Dy] = create_tv_operators(num_x, num_y);
    tv_kon = sum(abs(Dx*kon_params)) + sum(abs(Dy*kon_params));
    tv_koff = sum(abs(Dx*koff_params)) + sum(abs(Dy*koff_params));
    
    % Total cost
    misfit = norm(s_obs - s_fit);
    cost = misfit + lambda_kon*tv_kon + lambda_koff*tv_koff;
end

function [t, c_s, s] = simulate_3d_flow_model(nx, ny, nz, kon_grid,...
    koff_grid, smax_grid, velocity_profile, c0_assoc, c0_diss,...
    t_assoc, t_total, D_coeff, ru_to_m)
    
    % Initialize concentrations
    c_s = zeros(nx, ny, nz);
    c_s(1, :, :) = c0_assoc;  % Inlet at x=1
    s = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:)];

    % Time parameters
    tspan_assoc = linspace(0, t_assoc, 200);
    tspan_diss = linspace(t_assoc, t_total, 200);
    
    % Solve ODE
    options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz,...
        velocity_profile, kon_grid, koff_grid, smax_grid, c0_assoc,...
        D_coeff, ru_to_m), tspan_assoc, y0, options);
    
    % Reset for dissociation
    y_end_assoc = y_assoc(end,:)';
    c_s_end = reshape(y_end_assoc(1:nx*ny*nz), [nx, ny, nz]);
    s_end = reshape(y_end_assoc(nx*ny*nz+1:end), [nx, ny, nz]);  % Correct reshaping
    
    c_s_end(1, :, :) = c0_diss;  % Set inlet to 0
    y0_diss = [c_s_end(:); s_end(:)];  % Proper flattening of both components
    
    [t_diss, y_diss] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz,...
        velocity_profile, kon_grid, koff_grid, smax_grid, c0_diss,...
        D_coeff, ru_to_m), tspan_diss, y0_diss, options);
    
    % Combine results
    t = [t_assoc; t_diss(2:end)];
    y = [y_assoc; y_diss(2:end,:)];
    c_s = reshape(y(:,1:nx*ny*nz), [length(t), nx, ny, nz]);
    s = reshape(y(:,nx*ny*nz+1:end), [length(t), nx, ny, nz]);
end

function dydt = ode_system(t, y, nx, ny, nz, velocity_profile,...
    kon_grid, koff_grid, smax_grid, c0, D_coeff, ru_to_m)
    
    % Reshape state variables
    c_s = reshape(y(1:nx*ny*nz), [nx, ny, nz]);
    s = reshape(y(nx*ny*nz+1:end), [nx, ny, nz]);
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);
    
    % Diffusion terms
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:);
    
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2);
    d2c_dz2(:,:,1) = c_s(:,:,2) - 2*c_s(:,:,1) + c_s(:,:,1);
    d2c_dz2(:,:,end) = c_s(:,:,end-1) - 2*c_s(:,:,end) + c_s(:,:,end-1);
    
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2);
    
    % Advection
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)));
    
    % Adsorption kinetics
    dsdt = kon_grid .* c_s .* (smax_grid - s) - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c0;
    dcsdt(1,:,:) = 0;
    
    dydt = [dcsdt(:); dsdt(:)];
end