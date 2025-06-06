function heterogeneous_binding_analysis()
    % Main function for analyzing heterogeneous binding parameters
    clearvars; close all; clc;
    
    % ================ SIMULATION PARAMETERS ================
    gridN_x = 25;          % Reduced grid for faster testing
    gridN_y = 5;
    gridN_z = 3;
    ads_layer = 1;         % Adsorption layer (z=1)
    ads_x_range = [10,14]; % Adsorption region
    ads_y_range = [2,4];   % y-region
    c0_assoc = 3.3e-6;     % Association concentration
    c0_diss = 0;           % Dissociation concentration
    t_association = 800;   % Simulation times
    t_dissociation = 1600;
    D_coeff = 6e-3;        % Diffusion coefficient
    ru_to_m = 1e-6;        % RU conversion
    
    % Regularization parameters
    lambda_kon = 1e-2;
    lambda_koff = 1e-2;
    lambda_smax = 1e-2;
    lambda_mean = 1e2;     % Constraint strength
    
    % Flow parameters
    velocity_profile = create_velocity_profile(gridN_z, 8.3);

    % ================ SYNTHETIC DATA GENERATION ================
    fprintf('Generating synthetic data...\n');
    rng(1);  % For reproducibility
    
    % Create true parameter fields
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells)*length(ads_y_cells);
    [xx,yy] = meshgrid(ads_x_cells, ads_y_cells);
    
    % Smooth random parameter fields
    sz = [length(ads_x_cells), length(ads_y_cells)];
    x = linspace(-1,1,sz(1)); y = linspace(-1,1,sz(2));
    [xg,yg] = meshgrid(x,y); kernel = exp(-(xg.^2 + yg.^2));
    
    % Generate parameters with controlled variation
    noise_kon = randn(sz); noise_koff = randn(sz); noise_smax = randn(sz);
    
    % Base parameters (total smax = 1.0 RU)
    log_kon_base = log10(9.4e3);    % True mean kon
    log_koff_base = log10(0.0078);  % True mean koff
    smax_base = 1.0/num_ads_cells;  % Per-cell base capacity
    
    % Add correlated variations
    variation = 0.1; % 10% variation
    true_kon_log = log_kon_base + variation*imfilter(noise_kon, kernel, 'circular');
    true_koff_log = log_koff_base + variation*imfilter(noise_koff, kernel, 'circular');
    true_smax_grid = smax_base + variation*smax_base*imfilter(noise_smax, kernel, 'circular');
    
    % Convert to linear scale and flatten
    true_kon_grid = 10.^true_kon_log; true_koff_grid = 10.^true_koff_log;
    true_kon = true_kon_grid(:); true_koff = true_koff_grid(:); 
    true_smax = true_smax_grid(:);

    % Create 3D parameter grids
    [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        true_kon_grid, true_koff_grid, true_smax_grid);

    % Run forward model with updated parameters
    [t, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_grid,...
        koff_grid, smax_grid, velocity_profile, 0, 0, c0_assoc, c0_diss,...
        t_association, t_association + t_dissociation, ru_to_m, ads_layer,...
        1, [1, gridN_y], [1, gridN_z], gridN_x, [1, gridN_y], [1, gridN_z],...
        D_coeff, ads_x_range, ads_y_range);
    
    % Generate noisy observations
    s_obs = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    noise_level = 0.02;
    s_obs = s_obs + noise_level*max(s_obs)*randn(size(s_obs));

    % ================ TWO-STEP INVERSION ================
    % Step 1: Estimate homogeneous parameters
    fprintf('\n=== Step 1: Homogeneous parameter estimation ===\n');
    homog_params = estimate_homogeneous_parameters(t, s_obs, c0_assoc,...
        t_association, t_dissociation, ads_x_range, ads_y_range,...
        num_ads_cells, gridN_x, gridN_y, gridN_z, ads_layer,...
        velocity_profile, D_coeff, ru_to_m);
    
    % Step 2: Constrained heterogeneous inversion
    fprintf('\n=== Step 2: Heterogeneous inversion with constraints ===\n');
    [p_opt, animation_data] = solve_constrained_inversion(...
        homog_params, s_obs, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile,...
        c0_assoc, c0_diss, t_association, t_dissociation,...
        D_coeff, ru_to_m, lambda_kon, lambda_koff, lambda_smax, lambda_mean);

    % ================ VISUALIZATION ================
    generate_results_animation(animation_data, true_kon_grid,...
        true_koff_grid, true_smax_grid, t, s_obs);
end

%% Homogeneous parameter estimation (updated calls)
function homog_params = estimate_homogeneous_parameters(t, s_obs, c0_assoc,...
    t_assoc, t_diss, ads_x_range, ads_y_range, num_ads_cells,...
    gridN_x, gridN_y, gridN_z, ads_layer, velocity_profile, D_coeff, ru_to_m)
    
    % Parameter bounds [kon, koff, smax_total]
    lb = [1e3, 1e-5, 0.1];  
    ub = [1e5, 1e-1, 10];
    
    % Objective function with updated simulate_3d_flow_model call
    cost_func = @(p) homogeneous_cost(p, t, s_obs, c0_assoc, t_assoc, t_diss,...
        ads_x_range, ads_y_range, num_ads_cells, gridN_x, gridN_y, gridN_z,...
        ads_layer, velocity_profile, D_coeff, ru_to_m);
    
    % Optimization
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 10);
    homog_params = fmincon(cost_func, [9.4e3, 0.0078, 1.0], [], [], [], [], lb, ub, [], options);
    
    fprintf('Homogeneous parameters estimated:\n');
    fprintf('kon = %.2e M^{-1}s^{-1}\nkoff = %.2e s^{-1}\nsmax_total = %.2f RU\n',...
        homog_params(1), homog_params(2), homog_params(3));
end

function cost = homogeneous_cost(p, t, s_obs, c0_assoc, t_assoc, t_diss,...
    ads_x_range, ads_y_range, num_ads_cells, gridN_x, gridN_y, gridN_z,...
    ads_layer, velocity_profile, D_coeff, ru_to_m)
    
    % Create homogeneous parameter grids
    kon = p(1); 
    koff = p(2);
    smax_total = p(3);
    
    kon_grid = zeros(gridN_x, gridN_y, gridN_z);
    koff_grid = zeros(gridN_x, gridN_y, gridN_z);
    smax_grid = zeros(gridN_x, gridN_y, gridN_z);
    
    % Fill adsorption region
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_total/num_ads_cells;
    
    % Run simulation with updated parameters
    [~, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_grid,...
        koff_grid, smax_grid, velocity_profile, 0, 0, c0_assoc, 0,...
        t_assoc, t_assoc + t_diss, ru_to_m, ads_layer,...
        1, [1, gridN_y], [1, gridN_z], gridN_x, [1, gridN_y], [1, gridN_z],...
        D_coeff, ads_x_range, ads_y_range);
    
    % Calculate cost
    s_fit = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    cost = norm(s_obs - s_fit);
end

%% Constrained heterogeneous inversion (updated calls)
function [p_opt, animation_data] = solve_constrained_inversion(...
    homog_params, s_obs, gridN_x, gridN_y, gridN_z,...
    ads_x_range, ads_y_range, ads_layer, velocity_profile,...
    c0_assoc, c0_diss, t_assoc, t_total, D_coeff, ru_to_m,...
    lambda_kon, lambda_koff, lambda_smax, lambda_mean)
    
    % Parameter setup
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells)*length(ads_y_cells);
    
    % Initial guess (homogeneous params + noise)
    rng(2); % For reproducible initialization
    p0_kon = log10(homog_params(1)) + 0.5*randn(num_ads_cells,1);
    p0_koff = log10(homog_params(2)) + 0.5*randn(num_ads_cells,1);
    p0_smax = log10(homog_params(3)/num_ads_cells) + 0.5*randn(num_ads_cells,1);
    p0 = [p0_kon; p0_koff; p0_smax];
    
    % Bounds (log scale)
    lb_kon = log10(1e-2*ones(num_ads_cells,1));
    ub_kon = log10(1e5*ones(num_ads_cells,1));
    lb_koff = log10(1e-5*ones(num_ads_cells,1));
    ub_koff = log10(1e0*ones(num_ads_cells,1));
    lb_smax = log10(0.01*ones(num_ads_cells,1));
    ub_smax = log10(10*ones(num_ads_cells,1));
    lb = [lb_kon; lb_koff; lb_smax];
    ub = [ub_kon; ub_koff; ub_smax];
    
    % Initialize animation data
    global animation_data;
    animation_data = struct('iter', {}, 'p', {}, 's_fit', {});
    
    % Initial evaluation with updated simulate_3d_flow_model call
    [~, s_fit_initial] = constrained_cost(p0, s_obs, homog_params,...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff,...
        ru_to_m, lambda_kon, lambda_koff, lambda_smax, lambda_mean);
    
    animation_data(1).p = p0;
    animation_data(1).s_fit = s_fit_initial;
    animation_data(1).iter = 0;

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'iter',...
        'Algorithm', 'interior-point', 'MaxIterations', 100,...
        'UseParallel', true, 'OutputFcn', @outputfun,'MaxFunctionEvaluations',10000);
    
    % Run optimization
    p_opt = fmincon(@(p) constrained_cost(p, s_obs, homog_params,...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff,...
        ru_to_m, lambda_kon, lambda_koff, lambda_smax, lambda_mean),...
        p0, [], [], [], [], lb, ub, [], options);

    % Nested output function
    function stop = outputfun(p, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            [~, s_fit] = constrained_cost(p, s_obs, homog_params,...
                gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
                velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff,...
                ru_to_m, 0, 0, 0, 0); % No reg for visualization
            
            current_iter = optimValues.iteration + 1;
            animation_data(end+1).p = p;
            animation_data(end).s_fit = s_fit;
            animation_data(end).iter = current_iter;
        end
    end
end

function [cost, s_fit] = constrained_cost(p, s_obs, homog_params,...
    gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
    velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff,...
    ru_to_m, lambda_kon, lambda_koff, lambda_smax, lambda_mean)
    
    % Split parameters
    num_ads_cells = length(ads_x_range(1):ads_x_range(2)) * ...
                    length(ads_y_range(1):ads_y_range(2));
    kon_params = p(1:num_ads_cells);
    koff_params = p(num_ads_cells+1:2*num_ads_cells);
    smax_params = p(2*num_ads_cells+1:end);
    
    % Convert to linear scale
    kon_vals = 10.^kon_params;
    koff_vals = 10.^koff_params;
    smax_vals = 10.^smax_params;
    
    % Create parameter grids
    kon_grid = reshape(kon_vals, [length(ads_x_range(1):ads_x_range(2)),...
                                 length(ads_y_range(1):ads_y_range(2))]);
    koff_grid = reshape(koff_vals, size(kon_grid));
    smax_grid = reshape(smax_vals, size(kon_grid));
    
    [kon_3d, koff_3d, smax_3d] = generate_hetero_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        kon_grid, koff_grid, smax_grid);
    
    % Run simulation with updated parameters
    [~, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_3d,...
        koff_3d, smax_3d, velocity_profile, 0, 0, c0_assoc, c0_diss,...
        t_assoc, t_total, ru_to_m, ads_layer,...
        1, [1, gridN_y], [1, gridN_z], gridN_x, [1, gridN_y], [1, gridN_z],...
        D_coeff, ads_x_range, ads_y_range);
    
    % Calculate misfit
    s_fit = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    misfit = norm(s_obs - s_fit);
    
    % Regularization terms
    reg_kon = lambda_kon * norm(diff(kon_params))^2;
    reg_koff = lambda_koff * norm(diff(koff_params))^2;
    reg_smax = lambda_smax * norm(diff(smax_params))^2;
    
    % Constraint terms
    mean_kon = mean(kon_vals);
    mean_koff = mean(koff_vals);
    sum_smax = sum(smax_vals);
    
    target_kon = homog_params(1); 
    target_koff = homog_params(2);
    target_smax = homog_params(3);
    
    constraint_kon = lambda_mean * (mean_kon - target_kon)^2;
    constraint_koff = lambda_mean * (mean_koff - target_koff)^2;
    constraint_smax = lambda_mean * (sum_smax - target_smax)^2;
    
    % Total cost
    cost = misfit + reg_kon + reg_koff + reg_smax + ...
           constraint_kon + constraint_koff + constraint_smax;
end

%% Visualization function (unchanged)
function generate_results_animation(animation_data, true_kon_grid,...
    true_koff_grid, true_smax_grid, t, s_obs)
    
    fprintf('\nGenerating results animation...\n');
    video_filename = 'constrained_inversion_test.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 2;
    open(v);
    
    fig = figure('Position', [100 100 1200 1000]);
    
    for idx = 1:length(animation_data)
        current_p = animation_data(idx).p;
        current_s_fit = animation_data(idx).s_fit;
        iter_num = animation_data(idx).iter;
        
        num_ads_cells = numel(true_kon_grid);
        kon = 10.^current_p(1:num_ads_cells);
        koff = 10.^current_p(num_ads_cells+1:2*num_ads_cells);
        smax = 10.^current_p(2*num_ads_cells+1:end);
        
        kon_grid = reshape(kon, size(true_kon_grid));
        koff_grid = reshape(koff, size(true_koff_grid));
        smax_grid = reshape(smax, size(true_smax_grid));
        
        % Plotting
        clf(fig);
        colormap(gray);
        
        % True vs Recovered Parameters
        subplot(4,2,1);
        imagesc(log10(true_kon_grid'));
        colorbar; title('True log(k_{on})'); axis equal tight;
        
        subplot(4,2,2);
        imagesc(log10(kon_grid'));
        colorbar; title(sprintf('Recovered log(k_{on}) - Iter %d', iter_num));
        axis equal tight;
        
        subplot(4,2,3);
        imagesc(log10(true_koff_grid'));
        colorbar; title('True log(k_{off})'); axis equal tight;
        
        subplot(4,2,4);
        imagesc(log10(koff_grid'));
        colorbar; title(sprintf('Recovered log(k_{off}) - Iter %d', iter_num)); 
        axis equal tight;
        
        subplot(4,2,5);
        imagesc(true_smax_grid');
        colorbar; title('True s_{max}'); axis equal tight;
        
        subplot(4,2,6);
        imagesc(smax_grid');
        colorbar; title(sprintf('Recovered s_max - Iter %d', iter_num)); 
        axis equal tight;
        
        subplot(4,2,7:8);
        plot(t, s_obs, 'o', 'MarkerSize', 4, 'DisplayName', 'Observed');
        hold on;
        plot(t, current_s_fit, 'LineWidth', 2, 'DisplayName', 'Fit');
        xlabel('Time (s)'); ylabel('Response (RU)');
        legend('Location', 'best'); title('Sensorgram Comparison');
        hold off;
        
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    close(v);
    fprintf('Animation saved as %s\n', video_filename);
end

%% Helper Functions (unchanged)
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

%% Simulation Functions (Updated as per user's request)
function [t, c_s, s] = simulate_3d_flow_model(nx, ny, nz, kon_grid,...
    koff_grid, smax_grid, velocity_profile, k_flow_v, k_flow_d, c0_assoc, c0_diss,...
    t_assoc, t_total, ru_to_m, ads_layer, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Initialize state variables
    c_s = zeros(nx, ny, nz);
    c_s(inlet_x, inlet_y(1):inlet_y(2), inlet_z(1):inlet_z(2)) = c0_assoc;
    s = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:)];  

    % Time parameters
    tspan_assoc = linspace(0, t_assoc, 300);
    tspan_diss = linspace(t_assoc, t_total, 300);
    
    % Solve ODE with velocity_profile
    options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,velocity_profile,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_assoc,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_assoc, y0, options);
    
    % Reset for dissociation phase
    y_end_assoc = y_assoc(end,:)';
    c_s_end = reshape(y_end_assoc(1:nx*ny*nz), [nx, ny, nz]);
    s_end = reshape(y_end_assoc(nx*ny*nz+1:end), [nx, ny, nz]);
    c_s_end(inlet_x, inlet_y(1):inlet_y(2), inlet_z(1):inlet_z(2)) = c0_diss;
    y0_diss = [c_s_end(:); s_end(:)];

    [t_diss, y_diss] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,velocity_profile,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_diss,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_diss, y0_diss', options);
    
    % Combine results
    t = [t_assoc; t_diss(2:end)];
    y = [y_assoc; y_diss(2:end,:)];
    c_s = reshape(y(:,1:nx*ny*nz), [length(t), nx, ny, nz]);
    s = reshape(y(:,nx*ny*nz+1:end), [length(t), nx, ny, nz]);
end

function dydt = ode_system(t,y,nx,ny,nz,velocity_profile,~,~,kon_grid,...
    koff_grid,smax_grid,c0,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
    outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range)
    
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
    
    % Boundary conditions
    inlet_y_range = inlet_y(1):inlet_y(2);
    inlet_z_range = inlet_z(1):inlet_z(2);
    dcsdt(inlet_x, inlet_y_range, inlet_z_range) = 0;
    c_s(inlet_x, inlet_y_range, inlet_z_range) = c0;
    %dcsdt(outlet_x, :, :) = 0;
    
    % Adsorption kinetics
    ads_mask = zeros(nx, ny, nz);
    ads_mask(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = 1;
    active_kon = kon_grid .* ads_mask;
    active_koff = koff_grid .* ads_mask;
    active_smax = smax_grid .* ads_mask;
    %available_sites = max(active_smax - s, 0);
    dsdt = active_kon .* c_s .* (active_smax-s) - active_koff .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    dydt = [dcsdt(:); dsdt(:)];
end