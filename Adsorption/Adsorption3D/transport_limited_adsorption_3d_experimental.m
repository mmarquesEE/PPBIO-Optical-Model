function transport_limited_adsorption_3d_experimental
    % Main function for experimental 3D flow-driven adsorption
    clear all; close all; clc;
    % ==================== MOLECULE PARAMETERS ====================
    % ==================== SIMULATION PARAMETERS ====================
    gridN_x = 25;                  % Length dimension (11 mm)
    gridN_y = 5;                  % Depth dimension (1.7 mm)
    gridN_z = 3;                   % Height dimension (0.3 mm)
    roughness_scale = 0;           % No surface roughness
    c0 = 3.3e-6;                  % Inlet concentration (0.21 μM)
    c0_diss = 0;                   % Dissociation phase concentration
    k_flow_horizontal = 8.3;       % Horizontal flow rate (mm/s)
    k_flow_vertical = 0;           % No vertical flow
    k_flow_depth = 0;              % No depth flow
    ru_to_m = 1e-6;                % RU to molar conversion
    t_association = 800;           % Association phase duration
    t_dissociation = 1500;          % Total simulation time
    adsorption_z_layer = 1;        % Adsorption layer (z=1, bottom)
    D_coeff = 6e-3;                % Diffusion coefficient (mm²/s)

    % Localized inlet/outlet positions
    inlet_x = 1;                   % Inlet at first column
    inlet_y = [1, gridN_y];             % Inlet spans entire depth
    inlet_z = [1, gridN_z];              % Inlet spans all layers
    outlet_x = gridN_x;                 % Outlet at last column
    outlet_y = [1, gridN_y];            % Outlet spans entire depth
    outlet_z = [1, gridN_z];             % Outlet spans all layers
    ads_x_range = [round(gridN_x/2)-1,round(gridN_x/2)-1];          % Adsorption region x indices (center)
    ads_y_range = [round(gridN_y/2)-1,round(gridN_y/2)-1];         % Adsorption region y indices (1 mm width)
    
    h = gridN_z - 1; % Physical height (gridN_z layers span 0 to h)

    % Generate parabolic velocity profile
    z_indices = 0:(gridN_z - 1); % Physical z from 0 to h (MATLAB indices 1:gridN_z)
    velocity_profile = 4 * k_flow_horizontal * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1, 1, gridN_z]); % GPU array
    figure;
    plot(z_indices, squeeze(velocity_profile), 'LineWidth', 2);
    xlabel('Vertical Position (z)');
    ylabel('Horizontal Velocity (mm/s)');
    title('Parabolic Velocity Profile');
    grid on;
    % ==================== PARAMETER GENERATION ====================
    [kon_grid, koff_grid, smax_grid] = generate_3d_parameters(...
        gridN_x, gridN_y, gridN_z, roughness_scale);
    
    % ==================== SIMULATION EXECUTION ====================
    [t, c_s, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z,...
        kon_grid, koff_grid, smax_grid, velocity_profile, k_flow_vertical,...
        k_flow_depth, c0, c0_diss, t_association, t_dissociation, ru_to_m,...
        adsorption_z_layer, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, D_coeff, ads_x_range,ads_y_range);

    % ==================== VISUALIZATION ====================
    create_3d_flow_animation(t, c_s, gridN_x, gridN_y, gridN_z, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z);
    create_main_figures(t, s, kon_grid, koff_grid, smax_grid, gridN_x, gridN_y, ads_x_range,ads_y_range);
    create_flow_animation(t, c_s, s, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, gridN_z);
end

%% Parameter Generation
function [kon_grid, koff_grid, smax_grid] = generate_3d_parameters(...
    nx, ny, nz, roughness_scale)
    % Create 3D parameter grids with adsorption only at z=1 (bottom layer)
    [X,Y] = meshgrid(linspace(-1,1,nx), linspace(-1,1,ny));
    Z = exp(-(X.^2 + Y.^2)/0.3);
    
    % Generate roughness pattern
    roughness = imfilter(randn(nx,ny), Z, 'circular');
    roughness = roughness/max(abs(roughness(:)));
    
    % Base parameters
    koff = 0.0078; % Given koff (1/s)
    kon = 9.4e3;   % Given kon (m³/(s·mol)) = 9400 1/(M·s)
    smax = 1;       % RU/site

    % Initialize 3D grids
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);

    % Populate bottom layer (z=1)
    kon_grid(:,:,1) = kon * (1 + roughness_scale*roughness);
    koff_grid(:,:,1) = koff * (1 - 0.4*roughness_scale*roughness);
    smax_grid(:,:,1) = smax;
end

%% 3D Flow Simulation (Modified to use velocity_profile)
function [t, c_s, s] = simulate_3d_flow_model(nx, ny, nz, kon_grid,...
    koff_grid, smax_grid, velocity_profile, k_flow_v, k_flow_d, c0_assoc, c0_diss,...
    t_assoc, t_total, ru_to_m, ads_layer, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Initialize state variables
    c_s = zeros(nx, ny, nz);
    c_s(inlet_x, inlet_y(1):inlet_y(2), inlet_z(1):inlet_z(2)) = c0_assoc;
    s = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:)];  

    % Time parameters
    tspan_assoc = linspace(0, t_assoc, t_assoc/10);
    tspan_diss = linspace(t_assoc, t_total, t_assoc/10);
    
    % Solve ODE with velocity_profile
    options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,velocity_profile,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_assoc,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_assoc, y0, options);
    
    % ===== RESET INLET CONCENTRATION FOR DISSOCIATION =====
    % Extract final state from association
    y_end_assoc = y_assoc(end,:)';
    c_s_end = reshape(y_end_assoc(1:nx*ny*nz), [nx, ny, nz]);
    s_end = reshape(y_end_assoc(nx*ny*nz+1:end), [nx, ny, nz]);

    % Set inlet concentration to c0_diss (0)
    c_s_end(inlet_x, inlet_y(1):inlet_y(2), inlet_z(1):inlet_z(2)) = c0_diss;
    y0_diss = [c_s_end(:); s_end(:)];  % New initial condition

    [t_diss, y_diss] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,velocity_profile,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_diss,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_diss, y0_diss', options);
    
    % Combine results
    t = [t_assoc; t_diss];
    y = [y_assoc; y_diss];
    
    c_s = reshape(y(:,1:nx*ny*nz), [length(t), nx, ny, nz]);
    s = reshape(y(:,nx*ny*nz+1:end), [length(t), nx, ny, nz]);
end

%% ODE System (Updated)
function dydt = ode_system(t,y,nx,ny,nz,velocity_profile,~,~,kon_grid,...
    koff_grid,smax_grid,c0,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
    outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Reshape state variables
    c_s = reshape(y(1:nx*ny*nz), [nx, ny, nz]);
    s = reshape(y(nx*ny*nz+1:end), [nx, ny, nz]);
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);
    
    % ==================== DIFFUSION (x and z only) ====================
    % Second derivative in x
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:);
    
    % Second derivative in z (with no-flux BCs: ∂C/∂z=0 at z=0 and z=h)
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2);
    
    % No-flux at z=0 (mirror layer 2)
    d2c_dz2(:,:,1) = c_s(:,:,2) - 2*c_s(:,:,1) + c_s(:,:,1);
    % No-flux at z=h (mirror layer end-1)
    d2c_dz2(:,:,end) = c_s(:,:,end-1) - 2*c_s(:,:,end) + c_s(:,:,end-1);
    
    % Total diffusion
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2) + dcsdt;
    
    % ==================== ADVECTION (x-direction, parabolic profile) ====================
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)));
    
    % ==================== BOUNDARY CONDITIONS ====================
    % Inlet (Dirichlet: C(t,0,z) = c0)
    inlet_y_range = inlet_y(1):inlet_y(2);
    inlet_z_range = inlet_z(1):inlet_z(2);
    
    % Enforce fixed concentration at inlet
    dcsdt(inlet_x, inlet_y_range, inlet_z_range) = 0;  % No change over time
    c_s(inlet_x, inlet_y_range, inlet_z_range) = c0;    % Override concentration
    
    % Outlet (Neumann: ∂C/∂x = 0 at x=l)
    dcsdt(outlet_x, :, :) = 0;  % No change in concentration at outlet
    
    % ==================== ADSORPTION KINETICS ====================
    ads_mask = zeros(nx, ny, nz);
    ads_mask(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = 1;
    
    active_kon = kon_grid .* ads_mask;
    active_koff = koff_grid .* ads_mask;
    active_smax = smax_grid .* ads_mask;
    
    dsdt = active_kon .* c_s .* (active_smax - s) - active_koff .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    dydt = [dcsdt(:); dsdt(:)];
end


%% Updated create_flow_animation function
function create_flow_animation(t, c_s, s, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, gridN_z)
    filename = '3D_flow_animation_localized.gif';
    h = figure('Position', [100 100 1000 400], 'Color', 'w');
    
    max_conc = max(c_s(:));
    max_cov = max(s(:));
    
    for i = 1:length(t)
        figure(h); % Ensure this figure is the current one
        % Concentration slice (x-z plane at mid-depth)
        subplot(1,2,1)
        mid_y = round(size(c_s,3)/2);
        conc_slice = squeeze(c_s(i,:,mid_y,:))';
        imagesc(conc_slice)
        hold on;
        % Plot inlet and outlet
        plot(inlet_x, inlet_z, 'g.', 'MarkerSize', 20)
        plot(outlet_x, outlet_z, 'r.', 'MarkerSize', 20)
        hold off;
        title(sprintf('Conc. @ y=%d, t=%.1fs', mid_y, t(i)))
        xlabel('Length (x)'), ylabel('Height (z)')
        axis equal tight
        colorbar
        caxis([0 max_conc])
        
        % Surface coverage (x-y plane at z=0)
        subplot(1,2,2)
        surface_cov = squeeze(s(i,:,:,1));
        imagesc(surface_cov)
        title('Surface Coverage (z=0 plane)')
        ylabel('Length (x)'), xlabel('Depth (y)')
        axis equal tight
        colorbar
        caxis([0 max_cov])
        
        % Save frame to GIF
        frame = getframe(h);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if i == 1
            imwrite(imind, cm, filename, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end


%% Visualization Functions for 3D System
function create_main_figures(t, s, kon_grid, koff_grid, smax_grid, gridN_x, gridN_y, ads_x_range,ads_y_range)
    % Creates analysis figures for 3D flow model
    figure('Name','3D System Analysis','Position',[100 100 1200 400])
    
    % ==================== Parameter Correlations ====================
    subplot(1,3,1)
    % Extract surface parameters (z=0 layer)
    surface_kon = kon_grid(:,:,1);
    surface_koff = koff_grid(:,:,1);
    surface_smax = smax_grid(:,:,1);
    
    scatter3(surface_kon(:), surface_koff(:), surface_smax(:), 50, 'filled')
    title('Surface Parameter Correlations')
    xlabel('k_{on}'), ylabel('k_{off}'), zlabel('s_{max}')
    grid on
    
    % ==================== Total Coverage ====================
    subplot(1,3,2)
    % Sum over all x,y positions in surface layer (z=0)
    total_coverage = sum(s(:,:,:,1), [2 3]);  
    % Calculate total_smax for normalization
    num_ads_x = ads_x_range(2) - ads_x_range(1) + 1;
    num_ads_y = ads_y_range(2) - ads_y_range(1) + 1;
    total_ads_sites = num_ads_x * num_ads_y;
    total_smax = total_ads_sites * 1; % Each site has smax 100
    plot(t, total_coverage/total_smax, 'LineWidth',2)
    title('Total Surface Coverage')
    xlabel('Time (s)'), ylabel('Total RU (z=0 plane)')

    % Final surface coverage (x-y plane at z=1)
    subplot(1,3,3)
    [X,Y] = meshgrid(1:gridN_x, 1:gridN_y);
    % Proper dimensions: [nx, ny] data with [ny, nx] grid requires transpose
    surf(X', Y', squeeze(s(end,:,:,1)))
    title('Final Surface Coverage (z=0 plane)')
    xlabel('Length (x)'), ylabel('Depth (y)'), zlabel('Coverage (RU)')
end

%% Updated create_3d_flow_animation function
function create_3d_flow_animation(t, c_s, gridN_x, gridN_y, gridN_z, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z)
    filename = 'experimental_3d_flow.gif';
    h = figure('Position', [100 100 800 600], 'Color', 'w');
    
    max_conc = max(c_s(:));
    [X,Y,Z] = meshgrid(1:gridN_x, 1:gridN_y, 1:gridN_z);
    
    for i = 1:2:length(t) % Sample every 2nd frame
        figure(h); % Ensure this figure is the current one
        conc3d = squeeze(c_s(i,:,:,:));
        conc3d_permuted = permute(conc3d, [2 1 3]);
        
        % Slice visualization
        slice(X,Y,Z,conc3d_permuted, [1 gridN_x], round(gridN_y/2), [1 gridN_z])
        hold on;
        
        % Inlet/outlet markers
        [inX, inY, inZ] = meshgrid(inlet_x, inlet_y(1):inlet_y(2), inlet_z(1):inlet_z(2));
        scatter3(inX(:), inY(:), inZ(:), 50, 'g', 'filled')
        
        [outX, outY, outZ] = meshgrid(outlet_x, outlet_y(1):outlet_y(2), outlet_z(1):outlet_z(2));
        scatter3(outX(:), outY(:), outZ(:), 50, 'r', 'filled')
        
        hold off;
        shading interp
        colormap(jet)
        caxis([0 max_conc])
        title(sprintf('3D Concentration @ t=%.1fs', t(i)))
        xlabel('Length (x)'), ylabel('Depth (y)'), zlabel('Height (z)')
        view(45,30)
        axis tight
        
        % Save frame
        frame = getframe(h);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if i == 1
            imwrite(imind, cm, filename, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end