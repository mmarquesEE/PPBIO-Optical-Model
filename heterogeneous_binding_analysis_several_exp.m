function heterogeneous_binding_analysis()
    % Main function for analyzing heterogeneous binding parameters
    clearvars; close all; clc;
    if isempty(gcp('nocreate'))
        parpool; % Start parallel pool
    end
    
    % ================ SIMULATION PARAMETERS ================
    gridN_x = 25;          % Reduced grid for faster testing
    gridN_y = 5;
    gridN_z = 3;
    ads_layer = 1;         % Adsorption layer (z=1)
    ads_x_range = [10,14]; % Adsorption region
    ads_y_range = [2,4];   % y-region
    t_association = 800;   % Simulation times
    t_dissociation = 1200;
    D_coeff = 6e-3;        % Diffusion coefficient cm^2/s
    ru_to_m = 1e-6;        % RU conversion
    
    % Regularization parameters
    lambda_kon = 1e-2;
    lambda_koff = 1e-2;
    lambda_smax = 1e-2;
    lambda_mean = 10;     % Theorem constraint strength
    lambda_afm = 0.1;      % AFM prior strength
    
    % ================ MULTI-EXPERIMENT DESIGN ================
    experiments = [
        struct('c0_assoc', 3.3e-6, 'c0_diss', 0, 'vmax_factor', 1, 't_association', 800, 't_dissociation', 1200);
        struct('c0_assoc', 1e-5,   'c0_diss', 0, 'vmax_factor', 1, 't_association', 800, 't_dissociation', 1200);
        struct('c0_assoc', 0,      'c0_diss', 0, 'vmax_factor', 2, 't_association',800,   't_dissociation', 1200)
    ];
    num_experiments = length(experiments);
    
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
    
    % Base parameters
    log_kon_base = log10(9.4e3);
    log_koff_base = log10(0.0078);
    smax_base = 1.0/num_ads_cells;
    
    % Add correlated variations
    variation = 0.1;
    true_kon_log = log_kon_base + variation*imfilter(noise_kon, kernel, 'circular');
    true_koff_log = log_koff_base + variation*imfilter(noise_koff, kernel, 'circular');
    true_smax_grid = smax_base + variation*smax_base*imfilter(noise_smax, kernel, 'circular');
    
    % Convert to linear scale
    true_kon_grid = 10.^true_kon_log; 
    true_koff_grid = 10.^true_koff_log;
    true_smax = true_smax_grid(:);

    % Create 3D parameter grids
    [kon_grid, koff_grid, smax_grid] = generate_hetero_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        true_kon_grid, true_koff_grid, true_smax_grid);

    % ================ AFM PRIOR GENERATION ================
    afm_noise_level = 0.15;
    afm_kon = true_kon_grid .* (1 + afm_noise_level*randn(size(true_kon_grid)));
    afm_koff = true_koff_grid .* (1 + afm_noise_level*randn(size(true_koff_grid)));
    afm_smax = true_smax_grid .* (1 + afm_noise_level*randn(size(true_smax_grid)));
    
    % ================ GENERATE OBSERVATIONS ================
    s_obs_all = cell(num_experiments, 1);
    noise_level = 0.02;
    t_obs_all = cell(num_experiments, 1);
    parfor m = 1:3
        exp = experiments(m);
        velocity_profile = create_velocity_profile(gridN_z, 8.3 * exp.vmax_factor);
        [t, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_grid,...
            koff_grid, smax_grid, velocity_profile, exp.c0_assoc, exp.c0_diss,...
            exp.t_association, exp.t_dissociation, D_coeff, ru_to_m);
        
        obs = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
            ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
        s_obs_all{m} = obs + noise_level*max(obs)*randn(size(obs));
        t_obs_all{m} = t; % Store time vector
    end

    % ================ THEOREM-ENHANCED INVERSION ================
    fprintf('\n=== Two-step inversion with theorem constraints ===\n');
    
    % Step 1: Estimate homogeneous parameters using first experiment
    fprintf('\n=== Step 1: Homogeneous parameter estimation ===\n');
    homog_params = estimate_homogeneous_parameters(...
        experiments(1), s_obs_all{1}, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m);
    
    % Step 2: Theorem-constrained inversion
    fprintf('\n=== Step 2: Theorem-constrained inversion ===\n');
    [p_opt, animation_data] = solve_constrained_inversion(...
        homog_params, s_obs_all,t_obs_all, experiments, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,...
        lambda_kon, lambda_koff, lambda_smax, lambda_mean, lambda_afm,...
        afm_kon, afm_koff, afm_smax);

    % ================ VISUALIZATION ================
    generate_results_animation(animation_data, true_kon_grid,...
        true_koff_grid, true_smax_grid, s_obs_all, experiments);
end

%% Homogeneous parameter estimation
function homog_params = estimate_homogeneous_parameters(...
    exp, s_obs, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range,...
    ads_layer, D_coeff, ru_to_m)
    
    % Parameter bounds [kon, koff, smax_total]
    lb = [1e3, 1e-5, 0.1];  
    ub = [1e5, 1e-1, 2.0];
    
    velocity_profile = create_velocity_profile(gridN_z, 8.3 * exp.vmax_factor);
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Objective function
    cost_func = @(p) homogeneous_cost(p, s_obs, exp, gridN_x, gridN_y, gridN_z,...
        ads_x_range, ads_y_range, ads_layer, velocity_profile, D_coeff, ru_to_m);
    
    % Optimization
    options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel',true, 'MaxIterations',2);
    homog_params = fmincon(cost_func, [9.4e3, 0.0078, 1.0], [], [], [], [], lb, ub, [], options);
    
    fprintf('Homogeneous parameters estimated:\n');
    fprintf('kon = %.2e M^{-1}s^{-1}\nkoff = %.2e s^{-1}\nsmax_total = %.2f RU\n',...
        homog_params(1), homog_params(2), homog_params(3));
end

function cost = homogeneous_cost(p, s_obs, exp, gridN_x, gridN_y, gridN_z,...
    ads_x_range, ads_y_range, ads_layer, velocity_profile, D_coeff, ru_to_m)
    
    % Create homogeneous parameter grids
    kon = p(1); 
    koff = p(2);
    smax_total = p(3);
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    kon_grid = zeros(gridN_x, gridN_y, gridN_z);
    koff_grid = zeros(gridN_x, gridN_y, gridN_z);
    smax_grid = zeros(gridN_x, gridN_y, gridN_z);
    
    % Fill adsorption region
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_total/num_ads_cells;
    
    % Run simulation
    [~, ~, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z, kon_grid,...
        koff_grid, smax_grid, velocity_profile, exp.c0_assoc, exp.c0_diss,...
        exp.t_association, exp.t_dissociation, D_coeff, ru_to_m);
    
    % Calculate cost
    s_fit = squeeze(sum(s(:, ads_x_range(1):ads_x_range(2),...
        ads_y_range(1):ads_y_range(2), ads_layer), [2,3,4]));
    cost = norm(s_obs - s_fit);
end

%% Enhanced Constrained Inversion
function [p_opt, animation_data] = solve_constrained_inversion(...
    homog_params, s_obs_all,t_obs_all, experiments, gridN_x, gridN_y, gridN_z,...
    ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m,...
    lambda_kon, lambda_koff, lambda_smax, lambda_mean, lambda_afm,...
    afm_kon, afm_koff, afm_smax)
    
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
    
    % AFM prior in log space (flatten directly)
    afm_kon_flat = log10(afm_kon(:));
    afm_koff_flat = log10(afm_koff(:));
    afm_smax_flat = log10(afm_smax(:));
    afm_prior = [afm_kon_flat; afm_koff_flat; afm_smax_flat];
    
    % Initialize animation data
    global animation_data;
    animation_data = struct('iter', {}, 'p', {}, 's_fit', {}, 'params', {});
    
    % Optimization options
    options = optimoptions('fmincon', 'Display', 'iter',...
        'Algorithm', 'interior-point', 'MaxIterations', 15,...
        'UseParallel', true, 'OutputFcn', @outputfun, 'MaxFunctionEvaluations', 10000);
    
    % Run optimization
    p_opt = fmincon(@(p) theorem_constrained_cost(p, s_obs_all, t_obs_all,homog_params,...
        experiments, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        D_coeff, ru_to_m, lambda_kon, lambda_koff, lambda_smax, lambda_mean, lambda_afm, afm_prior),...
        p0, [], [], [], [], lb, ub, [], options);

    % Nested output function
    function stop = outputfun(p, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            current_iter = optimValues.iteration + 1;
            [~, s_fit, params] = theorem_constrained_cost(p, s_obs_all, t_obs_all,homog_params,...
                experiments, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
                D_coeff, ru_to_m, 0, 0, 0, 0, 0, afm_prior); % No reg for visualization
            
            animation_data(end+1).p = p;
            animation_data(end).s_fit = s_fit;
            animation_data(end).params = params;
            animation_data(end).iter = current_iter;
        end
    end
end

%% Precompute homogeneous kernels for all experiments
function K_homog_all = precompute_homogeneous_kernels(homog_params, experiments, ...
    gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m)
    
    num_experiments = length(experiments);
    K_homog_all = cell(num_experiments, 1);
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells)*length(ads_y_cells);
    
    for m = 1:num_experiments
        exp = experiments(m);
        velocity_profile = create_velocity_profile(gridN_z, 8.3 * exp.vmax_factor);
        
        % Create homogeneous grids
        kon_grid = repmat(homog_params(1), [gridN_x, gridN_y, gridN_z]);
        koff_grid = repmat(homog_params(2), [gridN_x, gridN_y, gridN_z]);
        smax_grid = zeros(gridN_x, gridN_y, gridN_z);
        smax_grid(ads_x_cells, ads_y_cells, ads_layer) = homog_params(3)/num_ads_cells;
        
        % Run simulation
        [~, ~, ~, K_heterog] = simulate_3d_flow_model(...
            gridN_x, gridN_y, gridN_z, kon_grid, koff_grid, smax_grid, ...
            velocity_profile, exp.c0_assoc, exp.c0_diss, ...
            exp.t_association, exp.t_dissociation, D_coeff, ru_to_m);
        
        % Extract adsorption region
        ads_cells = K_heterog(:, ads_x_cells, ads_y_cells, ads_layer);
        K_homog_all{m} = squeeze(sum(ads_cells, [2,3,4]));
    end
end

%% Theorem-Constrained Cost Function with multi-experiment support
function [cost, s_fit_all, params] = theorem_constrained_cost(p, s_obs_all, t_obs_all,homog_params, ...
    experiments, gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer, ...
    D_coeff, ru_to_m, lambda_kon, lambda_koff, lambda_smax, lambda_mean, lambda_afm, afm_prior)
    
    persistent K_homog_all
    num_experiments = length(experiments);
    ads_x_cells = ads_x_range(1):ads_x_range(2);
    ads_y_cells = ads_y_range(1):ads_y_range(2);
    num_ads_cells = length(ads_x_cells)*length(ads_y_cells);
    
    % Split parameters
    kon_params = p(1:num_ads_cells);
    koff_params = p(num_ads_cells+1:2*num_ads_cells);
    smax_params = p(2*num_ads_cells+1:end);
    
    kon_vals = 10.^kon_params;
    koff_vals = 10.^koff_params;
    smax_vals = 10.^smax_params;
    
    % Store parameters for visualization
    params.kon = reshape(kon_vals, [length(ads_x_cells), length(ads_y_cells)]);
    params.koff = reshape(koff_vals, size(params.kon));
    params.smax = reshape(smax_vals, size(params.kon));
    
    % Create 3D parameter grids
    kon_grid_2d = reshape(kon_vals, [length(ads_x_cells), length(ads_y_cells)]);
    koff_grid_2d = reshape(koff_vals, size(kon_grid_2d));
    smax_grid_2d = reshape(smax_vals, size(kon_grid_2d));
    
    [kon_3d, koff_3d, smax_3d] = generate_hetero_param_grids(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer,...
        kon_grid_2d, koff_grid_2d, smax_grid_2d);
    
    % Precompute homogeneous kernels once
    if isempty(K_homog_all)
        K_homog_all = precompute_homogeneous_kernels(homog_params, experiments, ...
            gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer, D_coeff, ru_to_m);
    end
    
    cost = 0;
    s_fit_all = cell(num_experiments, 1);
    alpha_vals = kon_vals .* smax_vals;
    
    % Process each experiment
    for m = 1:num_experiments
        exp = experiments(m);
        velocity_profile = create_velocity_profile(gridN_z, 8.3 * exp.vmax_factor);
        
        % Run simulation
        [t, ~, s, K_heterog] = simulate_3d_flow_model(...
            gridN_x, gridN_y, gridN_z, kon_3d, koff_3d, smax_3d, ...
            velocity_profile, exp.c0_assoc, exp.c0_diss, ...
            exp.t_association, exp.t_dissociation, D_coeff, ru_to_m);
        
        % Extract adsorption region
        s_fit_sim = squeeze(sum(s(:, ads_x_cells, ads_y_cells, ads_layer), [2,3,4]));
        s_fit = interp1(t, s_fit_sim, t_obs_all{m}, 'linear', 'extrap');
        s_fit_all{m} = s_fit;
        
        % Misfit term
        misfit = norm(s_obs_all{m} - s_fit);
        cost = cost + misfit;
        
        % Kernel constraint (Theorem 1)
        ads_cells = K_heterog(:, ads_x_cells, ads_y_cells, ads_layer);
        sum_alphaK_heterog = zeros(size(t));
        for i = 1:length(t)
            K_slice = squeeze(ads_cells(i,:,:));
            sum_alphaK_heterog(i) = sum(alpha_vals .* K_slice(:));
        end
        kernel_misfit = trapz(t, (sum_alphaK_heterog - K_homog_all{m}).^2);
        cost = cost + lambda_mean * kernel_misfit;
    end
    
    % Spatial regularization (Laplacian)
    reg_kon = lambda_kon * (norm(diff(params.kon, 1, 1), 'fro')^2 + norm(diff(params.kon, 1, 2), 'fro')^2);
    reg_koff = lambda_koff * (norm(diff(params.koff, 1, 1), 'fro')^2 + norm(diff(params.koff, 1, 2), 'fro')^2);
    reg_smax = lambda_smax * (norm(diff(params.smax, 1, 1), 'fro')^2 + norm(diff(params.smax, 1, 2), 'fro')^2);
    cost = cost + reg_kon + reg_koff + reg_smax;
    
    % AFM prior term
    current_params = [kon_params; koff_params; smax_params];
    prior_term = lambda_afm * norm(current_params - afm_prior)^2;
    cost = cost + prior_term;
end

%% Visualization function
function generate_results_animation(animation_data, true_kon_grid,...
    true_koff_grid, true_smax_grid, s_obs_all, experiments)
    
    fprintf('\nGenerating results animation...\n');
    video_filename = 'heterogeneous_inversion_otherother.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 1;
    open(v);
    
    fig = figure('Position', [100 100 1200 900], 'Color', 'w');
    num_experiments = length(experiments);
    
    for idx = 1:length(animation_data)
        current_p = animation_data(idx).p;
        s_fit_all = animation_data(idx).s_fit;
        params = animation_data(idx).params;
        iter_num = animation_data(idx).iter;
        
        % Extract true parameters for comparison
        ads_x_cells = 1:size(true_kon_grid,1);
        ads_y_cells = 1:size(true_kon_grid,2);
        
        % Plotting
        clf(fig);
        
        % True vs Recovered Parameters
        subplot(4,3,1);
        imagesc(true_kon_grid');
        colorbar; title('True k_{on}'); axis equal tight;
        set(gca, 'FontSize', 8);
        
        subplot(4,3,2);
        imagesc(true_koff_grid');
        colorbar; title('True k_{off}'); axis equal tight;
        set(gca, 'FontSize', 8);
        
        subplot(4,3,3);
        imagesc(true_smax_grid');
        colorbar; title('True s_{max}'); axis equal tight;
        set(gca, 'FontSize', 8);
        
        subplot(4,3,4);
        imagesc(params.kon');
        colorbar; title(sprintf('Recovered k_{on} - Iter %d', iter_num));
        axis equal tight;
        set(gca, 'FontSize', 8);
        
        subplot(4,3,5);
        imagesc(params.koff');
        colorbar; title(sprintf('Recovered k_{off} - Iter %d', iter_num)); 
        axis equal tight;
        set(gca, 'FontSize', 8);
        
        subplot(4,3,6);
        imagesc(params.smax');
        colorbar; title(sprintf('Recovered s_{max} - Iter %d', iter_num)); 
        axis equal tight;
        set(gca, 'FontSize', 8);
        
        % Sensorgram comparison
        colors = lines(num_experiments);
        max_time = 0;
        for m = 1:num_experiments
            max_time = max(max_time, length(s_obs_all{m}));
        end
        t_vec = 1:max_time;
        
        subplot(4,1,4);
        hold on;
        for m = 1:num_experiments
            obs = s_obs_all{m};
            fit = s_fit_all{m};
            t_obs = linspace(0, max_time, length(obs));
            t_fit = linspace(0, max_time, length(fit));
            
            plot(t_obs, obs, 'o', 'Color', colors(m,:), 'MarkerSize', 3, 'DisplayName', sprintf('Exp %d Obs', m));
            plot(t_fit, fit, '-', 'LineWidth', 1.5, 'Color', colors(m,:), 'DisplayName', sprintf('Exp %d Fit', m));
        end
        xlabel('Time (s)'); ylabel('Response (RU)');
        title('Multi-Experiment Sensorgram Fitting');
        legend('Location', 'bestoutside');
        set(gca, 'FontSize', 9);
        grid on;
        hold off;
        
        % Add iteration info
        annotation('textbox', [0.05, 0.95, 0.9, 0.05], 'String', ...
            sprintf('Iteration %d | Regularization: kon=%.1e, koff=%.1e, smax=%.1e, Theorem=%.1e, AFM=%.1e', ...
            iter_num, lambda_kon, lambda_koff, lambda_smax, lambda_mean, lambda_afm), ...
            'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');
        
        drawnow;
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

function [t, c_s, s, K] = simulate_3d_flow_model(nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff, ru_to_m)
    % Initialize concentrations and auxiliary variables
    c_s = zeros(nx, ny, nz);
    c_s(1, :, :) = c0_assoc;
    s = zeros(nx, ny, nz);
    Q = zeros(nx, ny, nz);
    R = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:); Q(:); R(:)];

    % Time parameters
    tspan_assoc = linspace(0, t_assoc, round(t_assoc*max(velocity_profile)));
    tspan_diss = linspace(t_assoc, t_total, round(t_assoc*max(velocity_profile)));

    
    % Solve ODE for association phase
    options = odeset('RelTol',1e-5, 'AbsTol',1e-7);
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
    y = y(:); % Ensure column vector    
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
