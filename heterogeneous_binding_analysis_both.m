function heterogeneous_binding_analysis()
    clearvars; close all; clc;
    
    % ================ SIMULATION PARAMETERS ================
    gridN_x = 25; ads_x_range = [10,15]; % Adsorption region
    gridN_y = 5;  ads_y_range = [2,4];
    gridN_z = 3;  ads_layer = 1;
    c0_assoc = 3.3e-6; t_association = 500;
    c0_diss = 0; t_dissociation = 700;
    D_coeff = 6e-3; ru_to_m = 1e-6;
    
    % ================ GROUND TRUTH PARAMETERS ================
    [true_kon_grid, true_koff_grid] = create_ground_truth_pattern(...
        ads_x_range, ads_y_range);
    
    % ================ SYNTHETIC DATA GENERATION ================
    velocity_profile = create_velocity_profile(gridN_z, 8.3);
    [kon_grid, koff_grid, smax_grid] = generate_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        true_kon_grid, true_koff_grid);
    
    [t, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_grid,...
        koff_grid, smax_grid, velocity_profile, c0_assoc, c0_diss,...
        t_association, t_dissociation, D_coeff, ru_to_m);
    
    s_obs = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    s_obs = s_obs + 0.02*max(s_obs)*randn(size(s_obs));

    % ================ INVERSION SETUP ================
    num_ads_x = length(ads_x_range(1):ads_x_range(2));
    num_ads_y = length(ads_y_range(1):ads_y_range(2));
    num_params = 2 * num_ads_x * num_ads_y; % kon + koff per cell
    
    % Initial guess: [log10(kon_values); log10(koff_values)]
    p0 = [log10(1e2*ones(num_ads_x*num_ads_y,1)); 
          log10(1e-2*ones(num_ads_x*num_ads_y,1))];
    
    % Bounds: kon (1e-1-1e4), koff (1e-4-1e1)
    lb = [log10(1e-1*ones(num_ads_x*num_ads_y,1)); 
          log10(1e-4*ones(num_ads_x*num_ads_y,1))];
    ub = [log10(1e4*ones(num_ads_x*num_ads_y,1)); 
          log10(1e1*ones(num_ads_x*num_ads_y,1))];
    
    % Regularization
    lambda_kon = 1e-2; 
    lambda_koff = 5e-3;
    
    options = optimoptions('fmincon', 'Display', 'iter',...
        'Algorithm', 'interior-point', 'MaxIterations', 100,...
        'UseParallel', true);
    
    % ================ SOLVE INVERSION ================
    p_opt = fmincon(@(p) joint_cost_function(p, s_obs, num_ads_x, num_ads_y,...
        lambda_kon, lambda_koff, gridN_x, gridN_y, gridN_z, ads_x_range,...
        ads_y_range, ads_layer, velocity_profile, c0_assoc, c0_diss,...
        t_association, t_dissociation, D_coeff, ru_to_m),...
        p0, [], [], [], [], lb, ub, [], options);
    
    % ================ RESULTS & VISUALIZATION ================
    [kon_est, koff_est] = unpack_parameters(p_opt, num_ads_x, num_ads_y);
    
    figure;
    subplot(2,2,1); imagesc(log10(true_kon_grid')); colorbar; title('True log(k_{on})');
    subplot(2,2,2); imagesc(log10(kon_est')); colorbar; title('Estimated log(k_{on})');
    subplot(2,2,3); imagesc(log10(true_koff_grid')); colorbar; title('True log(k_{off})');
    subplot(2,2,4); imagesc(log10(koff_est')); colorbar; title('Estimated log(k_{off})');
    
    % Plot sensorgram comparison
    [~, s_fit] = joint_cost_function(p_opt, s_obs, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, c0_assoc,...
        c0_diss, t_association, t_dissociation, D_coeff, ru_to_m, 0); % No regularization for final fit
    
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

function [kon_grid, koff_grid] = create_ground_truth_pattern(ads_x_range, ads_y_range)
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    [X,Y] = meshgrid(ads_x_cells, ads_y_cells);
    
    % kon pattern
    kon_grid = 1e2*ones(size(X));
    kon_grid(X > mean(ads_x_range) & Y > mean(ads_y_range)) = 1e3;
    kon_grid(X < mean(ads_x_range) & Y < mean(ads_y_range)) = 1e1;
    
    % koff pattern (anti-correlated with kon)
    koff_grid = 1e-2*ones(size(X));
    koff_grid(X > mean(ads_x_range) & Y > mean(ads_y_range)) = 1e-3;
    koff_grid(X < mean(ads_x_range) & Y < mean(ads_y_range)) = 1e-1;
end

%% Parameter Grid Generation (Updated)
function [kon_3d, koff_3d, smax_3d] = generate_param_grids(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_grid, koff_grid)
    
    kon_3d = zeros(nx, ny, nz);
    koff_3d = zeros(nx, ny, nz);
    smax_3d = zeros(nx, ny, nz);
    
    x_idx = ads_x_range(1):ads_x_range(2);
    y_idx = ads_y_range(1):ads_y_range(2);
    
    kon_3d(x_idx, y_idx, ads_layer) = kon_grid;
    koff_3d(x_idx, y_idx, ads_layer) = koff_grid;
    smax_3d(x_idx, y_idx, ads_layer) = 1;
end
%% Create TV operators for 2D grid
function [Dx, Dy] = create_tv_operators(nx, ny)
    % Horizontal differences
    Dx = spdiags([-ones(nx,1) ones(nx,1)], [0 1], nx, nx);
    Dx = kron(speye(ny), Dx);
    
    % Vertical differences
    Dy = spdiags([-ones(ny,1) ones(ny,1)], [0 1], ny, ny);
    Dy = kron(Dy, speye(nx));
end
%% Parameter Handling Utilities
function [kon_grid, koff_grid] = unpack_parameters(p, num_ads_x, num_ads_y)
    num_cells = num_ads_x * num_ads_y;
    kon_log = p(1:num_cells);
    koff_log = p(num_cells+1:end);
    
    kon_grid = reshape(10.^kon_log, num_ads_x, num_ads_y);
    koff_grid = reshape(10.^koff_log, num_ads_x, num_ads_y);
end
function [cost, s_fit] = joint_cost_function(p, s_obs, num_ads_x, num_ads_y,...
    lambda_kon, lambda_koff, nx, ny, nz, ads_x_range, ads_y_range, ads_layer,...
    velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff, ru_to_m)
    
    % Unpack parameters
    [kon_grid, koff_grid] = unpack_parameters(p, num_ads_x, num_ads_y);
    
    % Generate 3D parameter grids
    [kon_3d, koff_3d, smax_3d] = generate_param_grids(...
        nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_grid, koff_grid);
    
    % Run forward model
    [~, ~, s] = simulate_3d_flow_model(nx, ny, nz, kon_3d, koff_3d,...
        smax_3d, velocity_profile, c0_assoc, c0_diss, t_assoc,...
        t_total, D_coeff, ru_to_m);
    
    % Extract signal
    s_fit = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    
    % Calculate TV regularization
    [Dx, Dy] = create_tv_operators(num_ads_x, num_ads_y);
    kon_vec = kon_grid(:);
    koff_vec = koff_grid(:);
    
    tv_kon = sum(abs(Dx*kon_vec)) + sum(abs(Dy*kon_vec));
    tv_koff = sum(abs(Dx*koff_vec)) + sum(abs(Dy*koff_vec));
    
    % Total cost
    cost = norm(s_obs - s_fit) + lambda_kon*tv_kon + lambda_koff*tv_koff;
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