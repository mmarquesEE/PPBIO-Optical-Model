function transport_limited_adsorption_3d_experimental
    % Main function for experimental 3D flow-driven adsorption
    clear all; close all; clc;

    % ==================== SIMULATION PARAMETERS ====================
    gridN_x = 11;                  % Length dimension (11 mm)
    gridN_y = 17;                  % Depth dimension (1.7 mm)
    gridN_z = 3;                   % Height dimension (0.3 mm)
    roughness_scale = 0;           % No surface roughness
    c0 = 3.3e-6;                  % Inlet concentration (0.21 μM)
    c0_diss = 0;                   % Dissociation phase concentration
    k_flow_horizontal = 8.3;       % Horizontal flow rate (mm/s)
    k_flow_vertical = 0;           % No vertical flow
    k_flow_depth = 0;              % No depth flow
    ru_to_m = 1e-6;                % RU to molar conversion
    t_association = 500;           % Association phase duration
    t_dissociation = 800;          % Total simulation time
    adsorption_z_layer = 1;        % Adsorption layer (z=1, bottom)
    D_coeff = 6e-5;                % Diffusion coefficient (mm²/s)

    % Localized inlet/outlet positions
    inlet_x = 1;                   % Inlet at first column
    inlet_y = [1, 17];             % Inlet spans entire depth
    inlet_z = [1, 3];              % Inlet spans all layers
    outlet_x = 11;                 % Outlet at last column
    outlet_y = [1, 17];            % Outlet spans entire depth
    outlet_z = [1, 3];             % Outlet spans all layers
    ads_x_range = [6,6];          % Adsorption region x indices (center)
    ads_y_range = [9,9];         % Adsorption region y indices (1 mm width)

    % ==================== PARAMETER GENERATION ====================
    [kon_grid, koff_grid, smax_grid] = generate_3d_parameters(...
        gridN_x, gridN_y, gridN_z, roughness_scale);

    % ==================== SIMULATION EXECUTION ====================
    [t, c_s, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z,...
        kon_grid, koff_grid, smax_grid, k_flow_horizontal, k_flow_vertical,...
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

%% 3D Flow Simulation
function [t, c_s, s] = simulate_3d_flow_model(nx, ny, nz, kon_grid,...
    koff_grid, smax_grid, k_flow_h, k_flow_v, k_flow_d, c0_assoc, c0_diss,...
    t_assoc, t_total, ru_to_m, ads_layer, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Initialize state variables (x,y,z)
    c_s = zeros(nx, ny, nz);
    c_s(inlet_x, inlet_y(1):inlet_y(2), inlet_z(1):inlet_z(2)) = c0_assoc;
    s = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:)];  % Flatten for ODE solver
    
    % Time parameters
    tspan_assoc = linspace(0, t_assoc, t_assoc*10);
    tspan_diss = linspace(t_assoc, t_total, t_assoc*10);
    
    % Solve ODE
    options = odeset('RelTol',1e-5,'AbsTol',1e-8);

    % Association phase
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,k_flow_h,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_assoc,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_assoc, y0, options);
    
    % Dissociation phase
    [t_diss, y_diss] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,k_flow_h,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_diss,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_diss, y_assoc(end,:)', options);
    
    % Combine results
    t = [t_assoc; t_diss];
    y = [y_assoc; y_diss];
    
    % Reshape results
    c_s = reshape(y(:,1:nx*ny*nz), [length(t), nx, ny, nz]);
    s = reshape(y(:,nx*ny*nz+1:end), [length(t), nx, ny, nz]);
end

%% ODE System
function dydt = ode_system(t,y,nx,ny,nz,k_flow_h,k_flow_v,k_flow_d,kon_grid,...
    koff_grid,smax_grid,c0,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Reshape state variables
    c_s = reshape(y(1:nx*ny*nz), [nx, ny, nz]);
    s = reshape(y(nx*ny*nz+1:end), [nx, ny, nz]);
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);
    
    % ==================== 3D FLOW TRANSPORT ====================
    % Horizontal (x-direction)
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + k_flow_h*(c_s(1:end-1,:,:) - c_s(2:end,:,:));
    dcsdt(:,2:end,:) = dcsdt(:,2:end,:) + k_flow_h*(c_s(:,1:end-1,:) - c_s(:,1:end-1,:));

    
    % ==================== DIFFUSION ====================
    laplacian = 6 * del2(c_s);
    dcsdt = dcsdt + D_coeff * laplacian;
    
    % ==================== BOUNDARY CONDITIONS ====================
    % Inlet (x=1, all y, all z)
    inlet_y_range = inlet_y(1):inlet_y(2);
    inlet_z_range = inlet_z(1):inlet_z(2);
    dcsdt(inlet_x, inlet_y_range, inlet_z_range) = dcsdt(inlet_x, inlet_y_range, inlet_z_range) + ...
        k_flow_h * (c0 - c_s(inlet_x, inlet_y_range, inlet_z_range));
    
    % Outlet (x=end, all y, all z)
    outlet_y_range = outlet_y(1):outlet_y(2);
    outlet_z_range = outlet_z(1):outlet_z(2);
    dcsdt(outlet_x, outlet_y_range, outlet_z_range) = dcsdt(outlet_x, outlet_y_range, outlet_z_range) - ...
        k_flow_h * c_s(outlet_x, outlet_y_range, outlet_z_range);
    
    % ==================== ADSORPTION KINETICS ====================
    ads_mask = zeros(nx, ny, nz); % x, y, z order
    ads_mask(...
        ads_x_range(1):ads_x_range(2),... % x indices
        ads_y_range(1):ads_y_range(2),... % y indices
        ads_layer...                       % z index
    ) = 1;
    
    % Verify grid sizes match
    assert(isequal(size(kon_grid), size(ads_mask)),...
        'Parameter grid and mask dimension mismatch');
    
    active_kon = kon_grid .* ads_mask;
    active_koff = koff_grid .* ads_mask;
    active_smax = smax_grid .* ads_mask;
    
    dsdt = active_kon .* c_s .* (active_smax - s) - active_koff .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Flatten derivatives
    dydt = [dcsdt(:); dsdt(:)];
end

%% Visualization (Updated)
function create_flow_animation(t, c_s, s, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, gridN_z)
    filename = '3D_flow_animation_localized.gif';
    h = figure('Position', [100 100 1000 400], 'Color', 'w');
    
    max_conc = max(c_s(:));
    max_cov = max(s(:));
    
    for i = 1:length(t)
        % Concentration slice (x-z plane at mid-depth)
        subplot(1,2,1)
        mid_y = round(size(c_s,3)/2);
        conc_slice = squeeze(c_s(i,:,mid_y,:))';
        imagesc(conc_slice)
        hold on;
        % Plot inlet (x=1, z=1) and outlet (x=end, z=end)
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
    
    % ==================== Coverage Dynamics ====================
    subplot(1,4,1)
    hold on
    x_positions = round(linspace(1, gridN_x));
    % Plot coverage at different x positions at mid-depth
    for x = x_positions
        plot(t, s(:, :, x), 'LineWidth', 2)
    end
    title('Coverage Dynamics')
    xlabel('Time (s)'), ylabel('Coverage (RU)')
    
    % ==================== Parameter Correlations ====================
    subplot(1,4,2)
    % Extract surface parameters (z=0 layer)
    surface_kon = kon_grid(:,:,1);
    surface_koff = koff_grid(:,:,1);
    surface_smax = smax_grid(:,:,1);
    
    scatter3(surface_kon(:), surface_koff(:), surface_smax(:), 50, 'filled')
    title('Surface Parameter Correlations')
    xlabel('k_{on}'), ylabel('k_{off}'), zlabel('s_{max}')
    grid on
    
    % ==================== Total Coverage ====================
    subplot(1,4,3)
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
    subplot(1,4,4)
    [X,Y] = meshgrid(1:gridN_x, 1:gridN_y);
    % Proper dimensions: [nx, ny] data with [ny, nx] grid requires transpose
    surf(X', Y', squeeze(s(end,:,:,1)))
    title('Final Surface Coverage (z=0 plane)')
    xlabel('Length (x)'), ylabel('Depth (y)'), zlabel('Coverage (RU)')
end

function create_3d_flow_animation(t, c_s, gridN_x, gridN_y, gridN_z, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z)
    filename = 'experimental_3d_flow.gif';
    h = figure('Position', [100 100 800 600], 'Color', 'w');
    
    max_conc = max(c_s(:));
    [X,Y,Z] = meshgrid(1:gridN_x, 1:gridN_y, 1:gridN_z);
    
    for i = 1:2:length(t) % Sample every 2nd frame
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