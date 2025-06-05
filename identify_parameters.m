function validate_kernel_with_pulses()
    clearvars; close all; clc;
    % Shared parameters
    gridN_x = 10; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [5,5]; ads_y_range = [2,2];
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    % Add physical constants
    ru_to_m = 1e-6;  % Surface concentration conversion factor (critical!)
    D_coeff = 3.3e-6; % Diffusion coefficient
    
    % Base parameters
    kon_base = 9.4e3; 
    koff_base = 0.0078;
    base_max_velocity = 8.3; % mm/s
    
    % Ground truth parameters
    [kon_grid_gt, koff_grid_gt, smax_grid_gt] = create_ground_truth_heterogeneity(...
        gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    
    % Extract adsorption region parameters
    kon_ads_gt = kon_grid_gt(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    koff_ads_gt = koff_grid_gt(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    smax_ads_gt = smax_grid_gt(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    params_gt = [kon_ads_gt(:); koff_ads_gt(:); smax_ads_gt(:)];
    N_p = length(params_gt);
    
    % Example pulse parameters
    T1 = 200; T2 = 400; T3 = 800;
    c1 = 3.3e-6; c2 = 5e-5; c_diss = 0;
    
    % Experiment settings
    exp_settings = struct(...
        'pulse_times', { [T1, T2, T3], [T1, T2, T3], [T1, T2, T3] }, ...
        'pulse_concs', { [c1, 0, c2], [c1, 0, 2*c2], [2*c1, 0, c2] }, ...
        't_total', {1200, 1200, 1200}, ...
        'max_velocity', {base_max_velocity, 2*base_max_velocity, 0.5*base_max_velocity}, ...
        'c_diss', {c_diss, c_diss, c_diss} ...
    );
    
    % 1. Simulate unperturbed experiments
    s_obs_unpert = cell(1,length(exp_settings));
    t_unpert = cell(1,length(exp_settings));
    fprintf('Simulating unperturbed pulse experiments...\n');
    for m = 1:length(exp_settings)
        [t, s] = simulate_heterogeneous_with_pulses(...
            kon_grid_gt, koff_grid_gt, smax_grid_gt, gridN_x, gridN_y, gridN_z, exp_settings(m), ru_to_m, D_coeff);
        s_obs_unpert{m} = compute_composite_signal(s, ads_x_range, ads_y_range, ads_layer);
        t_unpert{m} = t;
    end

    % =====================================================================
    % ADDED PLOTTING SECTION
    % =====================================================================
    % Plot composite signal and concentration profiles
    figure('Name', 'Pulse Experiment Signals', 'Position', [100, 100, 1200, 800]);
    for m = 1:length(exp_settings)
        % Plot composite signal
        subplot(2, length(exp_settings), m);
        plot(t_unpert{m}, s_obs_unpert{m}, 'b-', 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Composite Signal');
        title(sprintf('Exp %d: Signal (Vel = %.1f mm/s)', m, exp_settings(m).max_velocity));
        grid on;
        
        % Add pulse transition markers
        pulse_times = exp_settings(m).pulse_times;
        ylims = get(gca, 'YLim');
        hold on;
        for i = 1:length(pulse_times)
            plot([pulse_times(i), pulse_times(i)], ylims, 'r--', 'LineWidth', 1);
        end
        hold off;
        
        % Plot concentration profile
        subplot(2, length(exp_settings), length(exp_settings) + m);
        c_in = zeros(size(t_unpert{m}));
        pulse_times = exp_settings(m).pulse_times;
        pulse_concs = exp_settings(m).pulse_concs;
        for i = 1:length(t_unpert{m})
            t_val = t_unpert{m}(i);
            if t_val < pulse_times(1)
                c_in(i) = pulse_concs(1);
            elseif t_val < pulse_times(2)
                c_in(i) = pulse_concs(2);
            elseif t_val < pulse_times(3)
                c_in(i) = pulse_concs(3);
            else
                c_in(i) = exp_settings(m).c_diss;
            end
        end
        plot(t_unpert{m}, c_in, 'r-', 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Concentration (M)');
        title(sprintf('Exp %d: Inlet Concentration', m));
        grid on;
        
        % Add pulse transition markers
        ylims = get(gca, 'YLim');
        hold on;
        for i = 1:length(pulse_times)
            plot([pulse_times(i), pulse_times(i)], ylims, 'k--', 'LineWidth', 1);
        end
        hold off;
    end

    % 2. Finite difference Jacobian
    dp_rel = 1e-5;  % Relative perturbation
    J = cell(1,length(exp_settings));  % Jacobian for each experiment
    
    fprintf('Computing Jacobian with pulses...\n');
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
            [t, s_pert] = simulate_heterogeneous_with_pulses(...
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
    title('Parameter Identifiability Analysis with Pulses');
end

%% Modified helper functions to include pulses
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

function [t, s] = simulate_heterogeneous_with_pulses(kon_grid, koff_grid, smax_grid, nx, ny, nz, exp_settings, ru_to_m, D_coeff)
    % Unpack settings
    pulse_times = exp_settings.pulse_times;
    pulse_concs = exp_settings.pulse_concs;
    t_total = exp_settings.t_total;
    max_velocity = exp_settings.max_velocity;
    c_diss = exp_settings.c_diss;
    
    velocity_profile = create_velocity_profile(nz, max_velocity);
    s0_grid = zeros(nx, ny, nz);  % Start with empty surface
    
    % Setup continuous time span
    tspan = linspace(0, t_total, 800);
    
    % Solve ODE with custom inlet concentration
    options = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [t, y] = ode15s(@(t,y) ode_pulse_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, pulse_times, pulse_concs, c_diss, ru_to_m, D_coeff), tspan, initial_state(nx, ny, nz));
    
    % Extract s from solution
    num_cells = nx * ny * nz;
    s = reshape(y(:,num_cells+1:2*num_cells), [length(t), nx, ny, nz]);
end

function y0 = initial_state(nx, ny, nz)
    % Initial state vector: c_s, s, Q, R all zeros except inlet
    num_cells = nx * ny * nz;
    c_s_init = zeros(nx, ny, nz);
    c_s_init(1,:,:) = 0;  % will be overwritten by inlet in ode
    s_init = zeros(nx, ny, nz);
    Q_init = zeros(nx, ny, nz);
    R_init = zeros(nx, ny, nz);
    y0 = [c_s_init(:); s_init(:); Q_init(:); R_init(:)];
end

function dydt = ode_pulse_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, pulse_times, pulse_concs, c_diss, ru_to_m, D_coeff)
    % Determine current inlet concentration
    if t < pulse_times(1)
        c_in = pulse_concs(1);
    elseif t < pulse_times(2)
        c_in = pulse_concs(2);
    elseif t < pulse_times(3)
        c_in = pulse_concs(3);
    else
        c_in = c_diss;
    end
    
    % Reshape state variables
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);
    Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);
    R = reshape(y(3*num_cells + 1:4*num_cells), [nx, ny, nz]);
    
    % Initialize derivatives
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
    dcsdt = dcsdt - (dsdt * ru_to_m);  % CRITICAL: Add conversion factor
    
    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c_in;
    dcsdt(1,:,:) = 0;
    
    % Compute dQ/dt and dR/dt
    dQdt = kon_grid .* c_s + koff_grid;
    dRdt = c_s .* exp(Q);

    % Combine all derivatives
    dydt = [dcsdt(:); dsdt(:); dQdt(:); dRdt(:)];
end

function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz-1);
    h = nz-1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1,1,nz]);
end