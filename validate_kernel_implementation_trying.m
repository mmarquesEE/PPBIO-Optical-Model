%% KERNEL VALIDATION SCRIPT
function validate_kernel_implementation()
    clearvars; close all; clc;
    % Shared parameters
    gridN_x = 10; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [5,5]; ads_y_range = [2,2];
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Base parameters for experiments
    kon_base = 9.4e3; 
    koff_base = 0.0078;
    kd_est = koff_base / kon_base;  % Estimated KD
    base_max_velocity = 8.3;        % Base flow velocity (mm/s)
    
    % Ground truth heterogeneous parameters
    [kon_grid_gt, koff_grid_gt, smax_grid_gt] = create_ground_truth_heterogeneity(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
    % Extract adsorption region parameters
    kon_ads_gt = kon_grid_gt(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    koff_ads_gt = koff_grid_gt(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    smax_ads_gt = smax_grid_gt(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    params_gt = [kon_ads_gt(:); koff_ads_gt(:); smax_ads_gt(:)];  % Flattened parameters 
    N_p = length(params_gt);  % Number of parameters 
    
    % Experiment settings with pulse experiment added
    T1 = 300; T2 = 600; T3 = 900; % Pulse times (seconds)
    c1 = 3.3e-6; c2 = 5e-6;      % Pulse concentrations
    
    % Define pulse concentration profile
    pulse_c_in = @(t) c1*(t < T1) + 0*(t >= T1 & t < T2) + c2*(t >= T2);
    
    % Experiment settings (struct array)
    exp_settings = struct(...
        'c_in_func', { ... 
            @(t) 3.3e-6*(t < 800), ...                 % Exp 1: Association only
            @(t) 1e-5*(t < 800) + 0*(t >= 800), ...     % Exp 2: Association + dissociation
            @(t) 5e-5*(t < 800) + 1e-5*(t >= 800), ...  % Exp 3: Association + different conc
            pulse_c_in ...                              % Exp 4: Pulse experiment
        }, ...
        't_total', {1200, 1200, 1200, T3}, ...          % Total simulation times
        'max_velocity', {base_max_velocity, ...          % Flow velocities
                         2*base_max_velocity, ...
                         0.5*base_max_velocity, ...
                         base_max_velocity} ...
    );
    
    % 1. Simulate unperturbed experiments
    s_obs_unpert = cell(1,length(exp_settings));
    t_unpert = cell(1,length(exp_settings));
    fprintf('Simulating unperturbed experiments...\n');
    for m = 1:length(exp_settings)
        [~, t, s, ~] = simulate_heterogeneous_kernel(...
            kon_grid_gt, koff_grid_gt, smax_grid_gt, gridN_x, gridN_y, gridN_z, exp_settings(m));
        s_obs_unpert{m} = compute_composite_signal(s, ads_x_range, ads_y_range, ads_layer);
        t_unpert{m} = t;
        
        % Plot each experiment
        figure;
        plot(t, s_obs_unpert{m}, 'b-', 'LineWidth', 2);
        xlabel('Time (s)'); ylabel('s_{obs}(t)');
        title(sprintf('Experiment %d: Composite Signal', m));
        grid on;
    end
    
    % 2. Finite difference Jacobian
    dp_rel = 1e-3;  % Relative perturbation
    J = cell(1,length(exp_settings));  % Jacobian for each experiment
    
    fprintf('Computing Jacobian...\n');
    for j = 1:N_p
        % Perturb j-th parameter
        params_pert = params_gt;
        params_pert(j) = params_gt(j) * (1 + dp_rel);
        dp_abs = params_gt(j) * dp_rel;
        
        % Reconstruct perturbed grids
        kon_ads_pert = reshape(params_pert(1:num_ads_cells), size(kon_ads_gt));
        koff_ads_pert = reshape(params_pert(num_ads_cells+1:2*num_ads_cells), size(koff_ads_gt));
        smax_ads_pert = reshape(params_pert(2*num_ads_cells+1:end), size(smax_ads_gt));
        
        kon_grid_pert = kon_grid_gt;
        kon_grid_pert(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_ads_pert;
        koff_grid_pert = koff_grid_gt;
        koff_grid_pert(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_ads_pert;
        smax_grid_pert = smax_grid_gt;
        smax_grid_pert(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_ads_pert;
        
        % Simulate perturbed system for all experiments
        for m = 1:length(exp_settings)
            [~, ~, s_pert, ~] = simulate_heterogeneous_kernel(...
                kon_grid_pert, koff_grid_pert, smax_grid_pert, gridN_x, gridN_y, gridN_z, exp_settings(m));
            s_obs_pert = compute_composite_signal(s_pert, ads_x_range, ads_y_range, ads_layer);
            
            % Compute Jacobian column
            if j == 1  % Initialize Jacobian on first parameter
                J{m} = zeros(length(t_unpert{m}), N_p);
            end
            J{m}(:, j) = (s_obs_pert - s_obs_unpert{m}) / dp_abs;
        end
    end
    
    % 3. Stack Jacobians and compute rank
    J_all = vertcat(J{:});  %  matrix
    fprintf('Computing SVD...\n');
    [U, S, V] = svd(J_all, 'econ');
    svals = diag(S);
    tolerance = 1e-6 * max(svals);
    rank_J = sum(svals > tolerance);
    
    fprintf('Rank of stacked Jacobian: %d/%d\n', rank_J, N_p);
    fprintf('Smallest singular value: %.2e\n', min(svals));
    
    % 4. Plot singular values
    figure;
    semilogy(svals, 'o-', 'LineWidth', 2);
    hold on;
    yline(tolerance, 'r--', 'Tolerance', 'LineWidth', 1.5);
    xlabel('Singular value index');
    ylabel('Singular value (log scale)');
    title(sprintf('Singular Values of Stacked Jacobian (Rank = %d/%d)', rank_J, N_p));
    grid on;
    legend('Singular values', 'Tolerance');

    % Identify nullspace parameters
    nullspace = V(:, rank_J+1:end);
    nullity = sum(abs(nullspace), 2);  % Importance of each parameter
    
    figure;
    bar(nullity);
    xlabel('Parameter index');
    ylabel('Null space contribution');
    title('Parameter Identifiability Analysis');
end

%% Helper Functions
function [kon_grid, koff_grid, smax_grid] = create_homogeneous_grids(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    % Create uniform parameter grids
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

function [kon_grid, koff_grid, smax_grid] = create_ground_truth_heterogeneity(nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    kon_base = 9.4e3;
    koff_base = 0.0078;
    smax_total = 1.0;
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Initialize grids
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    % Add 5% variation to adsorption region
    rng(42); % For reproducibility
    kon_vals = kon_base * (1 + 0.05*randn(ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1));
    koff_vals = koff_base * (1 + 0.05*randn(size(kon_vals)));
    smax_vals = (smax_total/num_ads_cells) * (1 + 0.05*randn(size(kon_vals)));
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_vals;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_vals;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_vals;
end

function s_obs = compute_composite_signal(s, ads_x_range, ads_y_range, ads_layer)
    ads_slab = s(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    s_obs = sum(ads_slab, [2,3,4]);  % Sum over x,y,z dimensions
end

function [K, t, s, Q] = simulate_heterogeneous_kernel(kon_grid, koff_grid, smax_grid, nx, ny, nz, exp_settings)
    % Unpack settings
    c_in_func = exp_settings.c_in_func;
    t_total = exp_settings.t_total;
    max_velocity = exp_settings.max_velocity;
    
    velocity_profile = create_velocity_profile(nz, max_velocity);
    s0_grid = zeros(nx, ny, nz);  % Start with empty surface
    
    [t, c_s, s, K, Q] = simulate_3d_flow_model(...
        nx, ny, nz, kon_grid, koff_grid, smax_grid,...
        velocity_profile, c_in_func, t_total, 3.3e-6, 1e-6, s0_grid);
end

function [t, c_s, s, K, Q] = simulate_3d_flow_model(nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, c_in_func, t_total, D_coeff, ru_to_m, s0_grid)
    % Initialize concentrations and auxiliary variables
    c_s = zeros(nx, ny, nz);
    c_s(1, :, :) = c_in_func(0);  % Initial inlet concentration
    s = s0_grid;  % Use provided initial surface state
    Q = zeros(nx, ny, nz);
    R = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:); Q(:); R(:)];

    % Time parameters - single phase simulation
    tspan = linspace(0, t_total, 400);
    
    % Solve ODE
    options = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [t, y] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c_in_func, D_coeff, ru_to_m), tspan, y0, options);
    
    % Extract variables
    num_cells = nx*ny*nz;
    c_s = reshape(y(:,1:num_cells), [length(t), nx, ny, nz]);
    s = reshape(y(:,num_cells+1:2*num_cells), [length(t), nx, ny, nz]);
    Q = reshape(y(:,2*num_cells+1:3*num_cells), [length(t), nx, ny, nz]);
    R = reshape(y(:,3*num_cells+1:end), [length(t), nx, ny, nz]);
    
    % Compute kernel K(x,y,t) = exp(-Q) .* R
    K = exp(-Q) .* R;
end

function dydt = ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c_in_func, D_coeff, ru_to_m)
    % Reshape state variables
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);
    Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);
    R = reshape(y(3*num_cells + 1:4*num_cells), [nx, ny, nz]);
    
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);

    % Diffusion terms
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = (c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:));
    
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2);
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2);
    
    % Advection
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)));
    
    % Adsorption kinetics with surface capacity constraint
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Apply time-dependent inlet boundary
    c_in = c_in_func(t);
    c_s(1,:,:) = c_in;
    dcsdt(1,:,:) = 0;
    
    % Compute dQ/dt and dR/dt
    dQdt = kon_grid .* c_s + koff_grid;
    dRdt = c_s .* exp(Q);
    
    % Combine all derivatives
    dydt = [dcsdt(:); dsdt(:); dQdt(:); dRdt(:)];
end