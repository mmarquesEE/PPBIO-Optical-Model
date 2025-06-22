function heterogeneous_binding_analysis()
    % Main function for analyzing heterogeneous binding parameters with experimental data (association only)
    clearvars; close all; clc;
    
    % ================ EXPERIMENTAL DATA LOADING ================
    fprintf('Loading experimental data...\n');
    
    % Define directory and load data
    dir_path = 'C:\REPOS\PPBIO-Optical-Model\Data_HCL\data_excel';
    file_list = dir(fullfile(dir_path, '*.xlsx'));
    num_files = 2;length(file_list);
    concentrations = [0.06, 0.08]; % Modify based on your actual concentrations
    
    % Preallocate structure for experimental data
    experimental_data = struct('t', cell(num_files,1), 's_obs', cell(num_files,1),...
                           'c0_assoc', cell(num_files,1));
    
    for i = 1:num_files
        data = readtable(fullfile(dir_path, file_list(i).name));
        experimental_data(i).t = data.Time / 1e3;       % Convert ms to seconds
        experimental_data(i).s_obs = data.RefractiveIndex / 1e3;  % Convert to RU
        experimental_data(i).c0_assoc = concentrations(i);
    end
    
    % ================ SIMULATION PARAMETERS ================
    gridN_x = 25;          % Reduced grid for faster testing
    gridN_y = 5;
    gridN_z = 3;
    ads_layer = 1;         % Adsorption layer (z=1)
    ads_x_range = [10,10]; % Central adsorption region
    ads_y_range = [2,2];   % Adsorption region in y
    D_coeff = 2.5e-3;        % Diffusion coefficient
    ru_to_m = 1e-6;        % RU conversion factor
    
    % Flow parameters
    velocity_profile = create_velocity_profile(gridN_z, 8.3);
    
    % ================ INVERSION SETUP ================
    fprintf('Setting up inverse problem...\n');
    
    % Adsorption region cells
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells) * length(ads_y_cells);
    
    % Initial parameter guess (log scale)
    p0_kon = log10(1e3 * ones(num_ads_cells,1)); % Reasonable guess for kon
    p0_koff = log10(1e-3 * ones(num_ads_cells,1)); % Reasonable guess for koff
    p0_smax = log10(1.2);
    p0 = [p0_kon; p0_koff; p0_smax];
    
    % Parameter bounds (log scale)
    lb_kon = log10(1e-1 * ones(num_ads_cells,1));  % Lower bound for kon
    ub_kon = log10(1e5 * ones(num_ads_cells,1));   % Upper bound for kon
    lb_koff = log10(1e-5 * ones(num_ads_cells,1)); % Lower bound for koff
    ub_koff = log10(1e-1 * ones(num_ads_cells,1)); % Upper bound for koff
    lb_smax = log10(0.1);
    ub_smax = log10(10);
    lb = [lb_kon; lb_koff; lb_smax];
    ub = [ub_kon; ub_koff; ub_smax];
    
    % Regularization parameters
    lambda_kon = 1e-3;
    lambda_koff = 1e-3;

    % ================ SOLVE INVERSION ================
    fprintf('Solving inverse problem...\n');
    options = optimoptions('fmincon', 'Display', 'iter',...
        'Algorithm', 'interior-point', 'MaxIterations', 100,...
        'UseParallel', true);
    
    p_opt = fmincon(@(p) cost_function(p, experimental_data, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, D_coeff, ru_to_m,...
        lambda_kon, lambda_koff), p0, [], [], [], [], lb, ub, [], options);
    
    % Split optimized parameters
    num_ads_cells = length(ads_x_cells) * length(ads_y_cells);
    recovered_kon = 10.^p_opt(1:num_ads_cells);
    recovered_koff = 10.^p_opt(num_ads_cells+1:end-1);
    recovered_smax = 10.^p_opt(end);
    recovered_kon_grid = reshape(recovered_kon, length(ads_x_cells), []);
    recovered_koff_grid = reshape(recovered_koff, length(ads_x_cells), []);

    % ================ RESULTS & VISUALIZATION ================
    fprintf('\nRecovery results:\n');
    fprintf('Recovered s_max: %.2f RU\n\n', recovered_smax);
    
    % Parameter maps
    figure;
    subplot(1,2,1);
    imagesc(log10(recovered_kon_grid'));
    colorbar; title('Recovered log(k_{on})'); axis equal tight;
    xlabel('X position'); ylabel('Y position');
    
    subplot(1,2,2);
    imagesc(log10(recovered_koff_grid'));
    colorbar; title('Recovered log(k_{off})'); axis equal tight;
    xlabel('X position'); ylabel('Y position');

    % Sensorgram comparison
    figure;
    hold on;
    [~, s_fit_all] = cost_function(p_opt, experimental_data, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, D_coeff, ru_to_m, 0, 0);
    
    colors = lines(length(experimental_data));
    for i = 1:length(experimental_data)
        data = experimental_data(i);
        plot(data.t, data.s_obs, 'o', 'Color', colors(i,:), 'MarkerSize', 4,...
            'DisplayName', sprintf('Observed (c=%.2f)', data.c0_assoc));
        plot(data.t, s_fit_all{i}, '-', 'Color', colors(i,:), 'LineWidth', 2,...
            'DisplayName', sprintf('Fit (c=%.2f)', data.c0_assoc));
    end
    xlabel('Time (s)'); ylabel('Response (RU)');
    legend('Location', 'best'); 
    title('Association Phase Fitting');
    grid on;
end

%% Helper Functions
function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz - 1);
    h = nz - 1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1, 1, nz]);
end

function [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_values, koff_values, smax_value)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    x_idx = ads_x_range(1):ads_x_range(2);
    y_idx = ads_y_range(1):ads_y_range(2);
    
    kon_grid(x_idx, y_idx, ads_layer) = kon_values;
    koff_grid(x_idx, y_idx, ads_layer) = koff_values;
    smax_grid(x_idx, y_idx, ads_layer) = smax_value;
end

function [Dx, Dy] = create_tv_operators(nx, ny)
    % Create Total Variation regularization operators
    Dx = spdiags([-ones(nx,1) ones(nx,1)], [0 1], nx, nx);
    Dx = kron(speye(ny), Dx);
    
    Dy = spdiags([-ones(ny,1) ones(ny,1)], [0 1], ny, ny);
    Dy = kron(Dy, speye(nx));
end

%% Cost Function with Regularization
function [cost, s_fit_all] = cost_function(p, experimental_data, nx, ny, nz, ads_x_range,...
    ads_y_range, ads_layer, velocity_profile, D_coeff, ru_to_m, lambda_kon, lambda_koff)
    
    % Split parameters
    num_ads_cells = length(ads_x_range(1):ads_x_range(2)) * length(ads_y_range(1):ads_y_range(2));
    kon_params = p(1:num_ads_cells);
    koff_params = p(num_ads_cells+1:end-1);
    smax_value = 10.^p(end);
    
    % Reshape to grids
    num_x = length(ads_x_range(1):ads_x_range(2));
    num_y = length(ads_y_range(1):ads_y_range(2));
    kon_grid = reshape(10.^kon_params, num_x, num_y);
    koff_grid = reshape(10.^koff_params, num_x, num_y);
    
    % Generate 3D parameter grids
    [kon_3d, koff_3d, smax_3d] = generate_hetero_param_grids(...
        nx, ny, nz, ads_x_range, ads_y_range, ads_layer,...
        kon_grid, koff_grid, smax_value);
    
    % Initialize total cost and storage
    cost = 0;
    s_fit_all = cell(length(experimental_data), 1);
    
    % Process each experimental dataset
    for i = 1:length(experimental_data)
        data = experimental_data(i);
        
        % Run forward model (association phase only)
        [t_sim, ~, s_sim] = simulate_3d_flow_model(nx, ny, nz, kon_3d, koff_3d,...
            smax_3d, velocity_profile, data.c0_assoc, D_coeff, ru_to_m, max(data.t));
        
        % Extract signal from adsorption region
        s_fit = squeeze(sum(s_sim(:, ads_x_range(1):ads_x_range(2),...
            ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
        
        % Interpolate to experimental time points
        s_fit_interp = interp1(t_sim, s_fit, data.t, 'pchip');
        
        % Calculate misfit
        misfit = norm(data.s_obs - s_fit_interp)^2;
        cost = cost + misfit;
        
        % Store for output
        s_fit_all{i} = s_fit_interp;
    end
    
    % Add regularization terms
    [Dx, Dy] = create_tv_operators(num_x, num_y);
    tv_kon = sum(sqrt((Dx*kon_params).^2 + (Dy*kon_params).^2));
    tv_koff = sum(sqrt((Dx*koff_params).^2 + (Dy*koff_params).^2));
    
    cost = cost + lambda_kon*tv_kon + lambda_koff*tv_koff;
end

%% Simulation Functions (Association Phase Only)
function [t, c_s, s] = simulate_3d_flow_model(nx, ny, nz, kon_grid, koff_grid, smax_grid,...
    velocity_profile, c0_assoc, D_coeff, ru_to_m, t_total)
    
    % Initialize concentrations
    c_s = zeros(nx, ny, nz);
    c_s(1, :, :) = c0_assoc;  % Inlet at x=1
    s = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:)];
    
    % Time parameters
    tspan = linspace(0, t_total, 200);
    
    % Solve ODE
    options = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [t, y] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz,...
        velocity_profile, kon_grid, koff_grid, smax_grid, c0_assoc,...
        D_coeff, ru_to_m), tspan, y0, options);
    
    % Reshape results
    c_s = reshape(y(:,1:nx*ny*nz), [length(t), nx, ny, nz]);
    s = reshape(y(:,nx*ny*nz+1:end), [length(t), nx, ny, nz]);
end

function dydt = ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid,...
    c0_assoc, D_coeff, ru_to_m)
    
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
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Boundary condition (constant concentration at inlet)
    c_s(1,:,:) = c0_assoc;
    dcsdt(1,:,:) = 0;
    
    dydt = [dcsdt(:); dsdt(:)];
end