%% Main Script for Testing the Equivalence Theorem in 3D Flow-Diffusion Model
clear all; close all; clc;

%% ================ SIMULATION PARAMETERS ================
% Grid dimensions
nx = 25;          % Reduced grid for faster testing
ny = 5;
nz = 3;
ads_layer = 1;    % Adsorption layer (z=1)

% Base parameters
c0_assoc = 3.3e-6;    % Association concentration (M)
c0_diss = 0;          % Dissociation concentration
t_association = 2000;  % Time units (e.g., seconds)
t_total = 4000;
D_coeff = 6e-3;       % Diffusion coefficient (cm²/s)
ru_to_m = 1e-6;       % RU conversion factor

% Flow velocity profile (parabolic in z-direction)
velocity_profile = create_velocity_profile(nz, 8.3);

%% ================ PARAMETER SETS ================
% Parameter Set 1 (p1: Random parameters with fixed total sums)
rng(42); % Seed for reproducibility

% Generate kon1 with total sum 9.4e3
kon1 = rand(nx, ny);
kon1 = kon1 / sum(kon1(:)) * 9.4e3; % Normalize sum to 9.4e3

% Generate koff1 with total sum 0.0078
koff1 = rand(nx, ny);
koff1 = koff1 / sum(koff1(:)) * 0.0078; % Normalize sum to 0.0078

% Generate smax1 with total sum 1.0
smax1 = rand(nx, ny);
smax1 = smax1 / sum(smax1(:)) * 1.0; % Normalize sum to 1.0

% Parameter Set 2 (p2: Permuted parameters from p1, maintaining sums)
perm_order = randperm(nx*ny); % Random permutation of site indices
kon2 = reshape(kon1(perm_order), nx, ny);
koff2 = reshape(koff1(perm_order), nx, ny);
smax2 = reshape(smax1(perm_order), nx, ny);


%% ================ GENERATE PARAMETER GRIDS ================
% For p1
[kon_grid_p1, koff_grid_p1, smax_grid_p1] = generate_param_grids(...
    nx, ny, nz, ads_layer, kon1, koff1, smax1);

% For p2
[kon_grid_p2, koff_grid_p2, smax_grid_p2] = generate_param_grids(...
    nx, ny, nz, ads_layer, kon2, koff2, smax2);

%% ================ RUN SIMULATIONS ================
% Simulate for p1
[t_p1, ~, s_p1] = simulate_3d_flow_model(...
    nx, ny, nz, kon_grid_p1, koff_grid_p1, smax_grid_p1, velocity_profile, ...
    c0_assoc, c0_diss, t_association, t_total, D_coeff, ru_to_m);

% Simulate for p2
[t_p2, ~, s_p2] = simulate_3d_flow_model(...
    nx, ny, nz, kon_grid_p2, koff_grid_p2, smax_grid_p2, velocity_profile, ...
    c0_assoc, c0_diss, t_association, t_total, D_coeff, ru_to_m);

%% ================ COMPUTE OBSERVED SIGNALS ================
sobs_p1 = squeeze(sum(s_p1, [2 3 4]));
sobs_p2 = squeeze(sum(s_p2, [2 3 4]));

%% ================ ANALYZE RESULTS ================
% Check equivalence between p1 and p2
max_diff_p1_p2 = max(abs(sobs_p1 - sobs_p2));
fprintf('Max difference (p1 vs p2): %.2e RU\n', max_diff_p1_p2);

%% ================ PLOT RESULTS ================
figure;
subplot(1,2,1);
plot(t_p1, sobs_p1, 'b', t_p2, sobs_p2, 'r--');
legend('p1', 'p2 (permuted)');
title('Equivalence Test: Permuted Parameters');
xlabel('Time (s)'); ylabel('s_{obs} (RU)');


%% Helper Functions
function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz - 1);
    h = nz - 1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1, 1, nz]);
end

function [kon_grid, koff_grid, smax_grid] = generate_param_grids(...
    nx, ny, nz, ads_layer, kon_params, koff_params, smax_params)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    % Assign parameters to each (x,y) in the adsorption layer
    kon_grid(:, :, ads_layer) = kon_params;
    koff_grid(:, :, ads_layer) = koff_params;
    smax_grid(:, :, ads_layer) = smax_params;
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
    tspan_assoc = linspace(0, t_assoc, round(t_assoc*max(velocity_profile)));
    tspan_diss = linspace(t_assoc, t_total, round(t_assoc*max(velocity_profile)));
    
    % Solve ODE
    options = odeset('RelTol',1e-9,'AbsTol',1e-9);
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
    
    dydt = [dcsdt(:); dsdt(:)];
end