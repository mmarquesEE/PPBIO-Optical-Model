function heterogeneous_binding_analysis()
    % Main function for analyzing heterogeneous binding parameters
    clearvars; close all; clc;
    
    % ================ SIMULATION PARAMETERS ================
    gridN_x = 25;          % Reduced grid for faster testing
    gridN_y = 5;
    gridN_z = 3;
    ads_layer = 1;         % Adsorption layer (z=1)
    ads_x_range = [10,14]; % Expanded adsorption region
    ads_y_range = [2,4];   % Expanded y-region
    c0_assoc = 3.3e-6;     % Association phase concentration
    c0_diss = 0;           % Dissociation phase concentration
    t_association = 1000;  % Reduced time for testing
    t_dissociation = 2000;
    D_coeff = 6e-3;        % Diffusion coefficient
    ru_to_m = 1e-6;        % RU conversion factor
    
    % Regularization parameters
    lambda_kon = 1e-2;
    lambda_koff = 1e-2;
    lambda_smax = 1e-2;
    
    % Flow parameters
    velocity_profile = create_velocity_profile(gridN_z, 8.3);
    
    % ================ HETEROGENEOUS PARAMETER SETUP ================
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells) * length(ads_y_cells);
    [xx,yy] = meshgrid(ads_x_cells, ads_y_cells);
    
        % ================ SYNTHETIC DATA GENERATION ================
    fprintf('Generating synthetic data...\n');
    rng(1);  % For reproducibility
    
    % Create smooth parameter fields with desired base values
    sz = [length(ads_x_cells), length(ads_y_cells)];
    x = linspace(-1,1,sz(1));
    y = linspace(-1,1,sz(2));
    [xg,yg] = meshgrid(x,y);
    kernel = exp(-(xg.^2 + yg.^2));
    
    % Generate correlated random fields with controlled variation
    noise_kon = randn(sz);
    noise_koff = randn(sz);
    noise_smax = randn(sz);
    
    % Set base values in log scale and add filtered noise
    log_kon_base = log10(9.4e3);    % ~4.0
    log_koff_base = log10(0.0078);  % ~-2.11
    smax_base = 1.0/length(ads_x_cells)*length(ads_y_cells);                % RU/site
    
    % Apply Gaussian smoothing with controlled standard deviation
    variation_kon = 0.01;  % Controls kon variation in log space
    variation_koff = 0.01; % Controls koff variation in log space
    variation_smax = 0.01; % Controls smax variation in linear space
    
    true_kon_log = log_kon_base + variation_kon * imfilter(noise_kon, kernel, 'circular');
    true_koff_log = log_koff_base + variation_koff * imfilter(noise_koff, kernel, 'circular');
    true_smax_grid = smax_base + variation_smax * imfilter(noise_smax, kernel, 'circular');
    
    % Convert to linear scale for kinetic parameters
    true_kon_grid = 10.^true_kon_log;
    true_koff_grid = 10.^true_koff_log;
    
    % Flatten for storage
    true_kon = true_kon_grid(:);
    true_koff = true_koff_grid(:);
    true_smax = true_smax_grid(:);
    
    % Create parameter grids
    [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        true_kon_grid, true_koff_grid, true_smax_grid);
    
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
    p0_kon = log10(true_kon_grid(:)) + 0.5*randn(num_ads_cells,1);
    p0_koff = log10(true_koff_grid(:)) + 0.5*randn(num_ads_cells,1);
    p0_smax = log10(true_smax_grid(:)) + 0.5*randn(num_ads_cells,1);
    p0 = [p0_kon; p0_koff; p0_smax];
    
    % Parameter bounds (log-scale)
    lb_kon = log10(1e-2 * ones(num_ads_cells,1));
    ub_kon = log10(1e5 * ones(num_ads_cells,1));
    lb_koff = log10(1e-5 * ones(num_ads_cells,1));
    ub_koff = log10(1e0 * ones(num_ads_cells,1));
    lb_smax = log10(0.1*ones(num_ads_cells,1));
    ub_smax = log10(10*ones(num_ads_cells,1));
    lb = [lb_kon; lb_koff; lb_smax];
    ub = [ub_kon; ub_koff; ub_smax];
    
    % ================ SOLVE INVERSION ================
    fprintf('Solving inverse problem...\n');
    
    % Initialize storage for animation
    global animation_data;
    animation_data = struct('iter', {}, 'p', {}, 's_fit', {});
    
    % Compute initial fit and store
    [~, s_fit_initial] = cost_function(p0, s_obs, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, c0_assoc,...
        c0_diss, t_association, t_dissociation, D_coeff, ru_to_m,...
        lambda_kon, lambda_koff, lambda_smax);
    animation_data(1).p = p0;
    animation_data(1).s_fit = s_fit_initial;
    animation_data(1).iter = 0;

    % Define output function
    function stop = outputfun(p, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            % Compute s_fit for current p without regularization
            [~, s_fit] = cost_function(p, s_obs, gridN_x, gridN_y, gridN_z,...
                ads_x_range, ads_y_range, ads_layer, velocity_profile, c0_assoc,...
                c0_diss, t_association, t_dissociation, D_coeff, ru_to_m,...
                0, 0, 0); % No regularization for visualization
            
            % Store current iteration data
            current_iter = optimValues.iteration + 1; % Iteration starts at 0
            animation_data(end+1).p = p;
            animation_data(end).s_fit = s_fit;
            animation_data(end).iter = current_iter;
        end
    end

    options = optimoptions('fmincon', 'Display', 'iter',...
        'Algorithm', 'interior-point', 'MaxIterations', 20,...
        'UseParallel', true, 'OutputFcn', @outputfun);
    
    p_opt = fmincon(@(p) cost_function(p, s_obs, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, c0_assoc,...
        c0_diss, t_association, t_dissociation, D_coeff, ru_to_m,...
        lambda_kon, lambda_koff, lambda_smax), p0, [], [], [], [], lb, ub, [], options);
    

    % ================ GENERATE ANIMATION ================
    fprintf('Generating animation...\n');
    
    % Create video writer
    video_filename = 'parameter_evolution.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 2; % Adjust frame rate as needed
    open(v);
    
    % Create figure for animation
    fig = figure('Position', [100 100 1200 1000]);
    
    for idx = 1:length(animation_data)
        % Extract current data
        current_p = animation_data(idx).p;
        current_s_fit = animation_data(idx).s_fit;
        iter_num = animation_data(idx).iter;
        
        % Split parameters
        num = num_ads_cells;
        recovered_kon = 10.^current_p(1:num);
        recovered_koff = 10.^current_p(num+1:2*num);
        recovered_smax = 10.^current_p(2*num+1:end);
        
        % Reshape to grids
        recovered_kon_grid = reshape(recovered_kon, size(true_kon_grid));
        recovered_koff_grid = reshape(recovered_koff, size(true_koff_grid));
        recovered_smax_grid = reshape(recovered_smax, size(true_smax_grid));
        
        % Plot parameter maps
        clf(fig); % Clear current figure
        colormap(gray);  % Apply grayscale to all subplots
        % True vs Recovered Parameters
        subplot(4,2,1);
        imagesc(log10(true_kon_grid'));
        colorbar; title('True log(k_{on})'); axis equal tight;
        
        subplot(4,2,2);
        imagesc(log10(recovered_kon_grid'));
        colorbar; title(sprintf('Recovered log(k_{on}) - Iter %d', iter_num));
        axis equal tight;
        
        subplot(4,2,3);
        imagesc(log10(true_koff_grid'));
        colorbar; title('True log(k_{off})'); axis equal tight;
        
        subplot(4,2,4);
        imagesc(log10(recovered_koff_grid'));
        colorbar; title(sprintf('Recovered log(k_{off}) - Iter %d', iter_num)); 
        axis equal tight;
        
        subplot(4,2,5);
        imagesc(true_smax_grid');
        colorbar; title('True s_{max}'); axis equal tight;
        
        subplot(4,2,6);
        imagesc(recovered_smax_grid');
        colorbar; title(sprintf('Recovered s_max - Iter %d', iter_num)); 
        axis equal tight;
        
        % Sensorgram comparison in the fourth row spanning two columns
        subplot(4,2,7:8);
        plot(t, s_obs, 'o', 'MarkerSize', 4, 'DisplayName', 'Observed');
        hold on;
        plot(t, current_s_fit, 'LineWidth', 2, 'DisplayName', 'Fit');
        xlabel('Time (s)'); ylabel('Response (RU)');
        legend('Location', 'best'); title('Sensorgram Comparison');
        hold off;
        
        drawnow;
        
        % Capture frame
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    close(v);
    fprintf('Animation saved as %s\n', video_filename);
end

%% Helper Functions
function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz - 1);
    h = nz - 1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1, 1, nz]);
end

function [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_values, koff_values, smax_values)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    x_idx = ads_x_range(1):ads_x_range(2);
    y_idx = ads_y_range(1):ads_y_range(2);
    
    kon_grid(x_idx, y_idx, ads_layer) = kon_values;
    koff_grid(x_idx, y_idx, ads_layer) = koff_values;
    smax_grid(x_idx, y_idx, ads_layer) = smax_values;
end

%% Cost Function with Regularization
function [cost, s_fit] = cost_function(p, s_obs, nx, ny, nz, ads_x_range,...
    ads_y_range, ads_layer, velocity_profile, c0_assoc, c0_diss,...
    t_assoc, t_total, D_coeff, ru_to_m, lambda_kon, lambda_koff, lambda_smax)
    
    % Split parameters
    num_ads_cells = length(ads_x_range(1):ads_x_range(2)) * length(ads_y_range(1):ads_y_range(2));
    kon_params = p(1:num_ads_cells);
    koff_params = p(num_ads_cells+1:2*num_ads_cells);
    smax_params = p(2*num_ads_cells+1:end);
    
    % Reshape to grids
    num_x = length(ads_x_range(1):ads_x_range(2));
    num_y = length(ads_y_range(1):ads_y_range(2));
    kon_grid = reshape(10.^kon_params, num_x, num_y);
    koff_grid = reshape(10.^koff_params, num_x, num_y);
    smax_grid = reshape(10.^smax_params, num_x, num_y);
    
    % Generate 3D parameter grids
    [kon_3d, koff_3d, smax_3d] = generate_hetero_param_grids(...
        nx, ny, nz, ads_x_range, ads_y_range, ads_layer,...
        kon_grid, koff_grid, smax_grid);
    
    % Run forward model
    [~, ~, s] = simulate_3d_flow_model(nx, ny, nz, kon_3d, koff_3d,...
        smax_3d, velocity_profile, c0_assoc, c0_diss, t_assoc,...
        t_total, D_coeff, ru_to_m);
    
    % Extract signal
    s_fit = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    
    % Total cost
    misfit = norm(s_obs - s_fit);
    cost = misfit;
end

%% Simulation Functions (Unchanged from original)
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
    s_end = reshape(y_end_assoc(nx*ny*nz+1:end), [nx, ny, nz]);
    
    c_s_end(1, :, :) = c0_diss;
    y0_diss = [c_s_end(:); s_end(:)];
    
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
    
    % Adsorption kinetics with surface capacity constraint
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c0;
    dcsdt(1,:,:) = 0;
    
    dydt = [dcsdt(:); dsdt(:)];
end