function transport_limited_adsorption_3d_localized
    % Main function for 3D flow-driven adsorption with localized inlet/outlet
    clear all; close all; clc;

    % ==================== SIMULATION PARAMETERS ====================
    gridN_x = 10;                  % Length dimension (columns)
    gridN_y = 10;                  % Depth dimension (slices)
    gridN_z = 10;                  % Height dimension (rows)
    roughness_scale = 1;          % Surface roughness intensity (0-1)
    c0 = 200e-9;                  % Inlet concentration (M)
    c0_diss = 0;                  % Dissociation phase concentration
    k_flow_horizontal = 10;      % Horizontal flow rate (columns/s)
    k_flow_vertical = 10;        % Vertical flow rate (rows/s)
    k_flow_depth = 10;           % Depth flow rate (slices/s) 
    ru_to_m = 1e-6;               % RU to molar conversion
    t_association = 100;          % Association phase duration
    t_dissociation = 200;         % Total simulation time
    adsorption_z_layer = 1;       % Adsorption layer (z=1, bottom)
    D_coeff = 1e-3;               % Diffusion coefficient (cell²/s) [ADDED]

    % Localized inlet/outlet positions [NEW]
    inlet_x = 1;                  % Inlet at first column
    inlet_y = [1,1];             % Inlet spans depth y=3 to y=5
    inlet_z = gridN_z;                  % Inlet at bottom layer
    outlet_x = gridN_x;           % Outlet at last column
    outlet_y = [1,1];            % Outlet spans depth y=3 to y=5
    outlet_z = gridN_z;           % Outlet at top layer
    ads_x_range = [3, 8];          % Adsorption region x indices (columns)
    ads_y_range = [3, 8];          % Adsorption region y indices (depth)
    % ==================== PARAMETER GENERATION ====================
    [kon_grid, koff_grid, smax_grid] = generate_3d_parameters(...
    gridN_x, gridN_y, gridN_z, roughness_scale); % Added ads ranges

    % ==================== SIMULATION EXECUTION ====================
    [t, c_s, s] = simulate_3d_flow_model(gridN_x, gridN_y, gridN_z,...
        kon_grid, koff_grid, smax_grid, k_flow_horizontal, k_flow_vertical,...
        k_flow_depth, c0, c0_diss, t_association, t_dissociation, ru_to_m,...
        adsorption_z_layer, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, D_coeff, ads_y_range,ads_x_range);

    % ==================== VISUALIZATION ====================
    create_3d_flow_animation(t, c_s, gridN_x, gridN_y, gridN_z, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z); % New 3D animation
    create_main_figures(t, s, kon_grid, koff_grid, smax_grid, gridN_x, gridN_y, ads_y_range,ads_x_range);
    create_flow_animation(t, c_s, s, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, gridN_z);

end

%% Parameter Generation
function [kon_grid, koff_grid, smax_grid] = generate_3d_parameters(...
    nx, ny, nz, roughness_scale) % Modified input
    % Create 3D parameter grids with adsorption only at z=0 (bottom layer)
    [X,Y] = meshgrid(linspace(-1,1,nx), linspace(-1,1,ny));
    Z = exp(-(X.^2 + Y.^2)/0.3);
    
    % Generate roughness pattern for surface (z=0)
    roughness = imfilter(randn(nx,ny), Z, 'circular');
    roughness = roughness/max(abs(roughness(:)));
    
    % Base parameters
    koff = 1e-3; 
    KD = 1e-9;
    % Initialize 3D grids
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    % Only populate bottom layer (z=0)
    kon_grid(:,:,1) = (koff/KD) * (1 + roughness_scale*roughness);
    koff_grid(:,:,1) = koff * (1 - 0.4*roughness_scale*roughness);
    smax_grid(:,:,1) = 1;  % RU/site
end

%% 3D Flow Simulation (Updated)
function [t, c_s, s] = simulate_3d_flow_model(nx, ny, nz, kon_grid,...
    koff_grid, smax_grid, k_flow_h, k_flow_v, k_flow_d, c0_assoc, c0_diss,...
    t_assoc, t_total, ru_to_m, ads_layer, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Initialize state variables (x,y,z)
    c_s = zeros(nx, ny, nz);
    c_s(inlet_x, inlet_y, inlet_z) = c0_assoc;
    s = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:)];  % Flatten for ODE solver
    % Time parameters
    tspan_assoc = linspace(0, t_assoc,t_assoc);
    tspan_diss = linspace(t_assoc, t_total, t_assoc);
    
    % Solve ODE
    options = odeset('RelTol',1e-5,'AbsTol',1e-8);

    % Association phase
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,k_flow_h,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_assoc,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_assoc, y0, options); % [MODIFIED]
    
    % Dissociation phase
    [t_diss, y_diss] = ode15s(@(t,y) ode_system(t,y,nx,ny,nz,k_flow_h,k_flow_v,k_flow_d,...
        kon_grid,koff_grid,smax_grid,c0_diss,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,...
        outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range), tspan_diss, y_assoc(end,:)', options); % [MODIFIED]
    
    % Combine results
    t = [t_assoc; t_diss];
    y = [y_assoc; y_diss];
    
    % Reshape results
    c_s = reshape(y(:,1:nx*ny*nz), [length(t), nx, ny, nz]);
    s = reshape(y(:,nx*ny*nz+1:end), [length(t), nx, ny, nz]);
end

%% ODE System (Updated for 3D flow)
function dydt = ode_system(t,y,nx,ny,nz,k_flow_h,k_flow_v,k_flow_d,kon_grid,...
    koff_grid,smax_grid,c0,ru_to_m,ads_layer,inlet_x,inlet_y,inlet_z,outlet_x,outlet_y,outlet_z, D_coeff, ads_x_range, ads_y_range)
    
    % Reshape state variables
    c_s = reshape(y(1:nx*ny*nz), [nx, ny, nz]);
    s = reshape(y(nx*ny*nz+1:end), [nx, ny, nz]);
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);
    
    % ==================== 3D FLOW TRANSPORT ====================
    % Horizontal (x-direction)
    for x = 2:nx
        dcsdt(x,:,:) = dcsdt(x,:,:) + k_flow_h*(c_s(x-1,:,:) - c_s(x,:,:));
    end
    
    % Depth (y-direction) [NEW]
    for y = 2:ny
        dcsdt(:,y,:) = dcsdt(:,y,:) + k_flow_d*(c_s(:,y-1,:) - c_s(:,y,:));
    end
    
    % Vertical (z-direction) [CORRECTED DOWNWARD FLOW]
    for z = 1:nz-1
        dcsdt(:,:,z) = dcsdt(:,:,z) + k_flow_v*(c_s(:,:,z+1) - c_s(:,:,z));
    end
    % ==================== DIFFUSION ====================
    % Compute Laplacian for diffusion with Neumann BCs
    laplacian = 6 * del2(c_s);  % Adjust scaling for 3D
    dcsdt = dcsdt + D_coeff * laplacian; % [ADDED]
   % ==================== BOUNDARY CONDITIONS ====================
    % Localized inlet (x=1, y=1:8, z=top)
    inlet_y_range = inlet_y(1):inlet_y(2);
    dcsdt(inlet_x, inlet_y_range, inlet_z) = dcsdt(inlet_x, inlet_y_range, inlet_z) + ...
        k_flow_h *(c0 - c_s(inlet_x, inlet_y_range, inlet_z));
    
    % Localized outlet (x=end, y=1:8, z=1) [CORRECTED TO HORIZONTAL FLOW]
    outlet_y_range = outlet_y(1):outlet_y(2);
    dcsdt(outlet_x, outlet_y_range, outlet_z) = dcsdt(outlet_x, outlet_y_range, outlet_z) - ...
        k_flow_h  * c_s(outlet_x, outlet_y_range, outlet_z);
    
    % ==================== ADSORPTION KINETICS ====================
    ads_mask = zeros(nx,ny,nz);
    ads_mask(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2),ads_layer) = 1;  % Only at z=ads_layer
    active_kon = kon_grid .* ads_mask;
    active_koff = koff_grid .* ads_mask;
    active_smax = smax_grid .* ads_mask;
    
    dsdt = active_kon .* c_s .* (active_smax - s) - active_koff .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);  % Convert RU/s to concentration
    
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
        xlabel('Length (x)'), ylabel('Depth (y)')
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
function create_main_figures(t, s, kon_grid, koff_grid, smax_grid, gridN_x, gridN_y, ads_y_range,ads_x_range)
    % Creates analysis figures for 3D flow model
    figure('Name','3D System Analysis','Position',[100 100 1200 400])
    
    % ==================== Coverage Dynamics ====================
    subplot(1,4,1)
    hold on
    x_positions = round(linspace(1, gridN_x));
    center_size = round(gridN_y * 1);
    start_idx = floor((gridN_y - center_size)/2) + 1;
    end_idx = start_idx + center_size - 1;
    rows_to_plot = round(linspace(start_idx,end_idx));
    % Plot coverage at different x positions at mid-depth
    for x = x_positions
        plot(t, mean(s(:, rows_to_plot, x), 2), 'LineWidth', 2)
    end
    title('Coverage Dynamics by Column')
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

%% Corrected 3D Fluid Flow Animation Function
function create_3d_flow_animation(t, c_s, gridN_x, gridN_y, gridN_z, inlet_x, inlet_y, inlet_z, outlet_x, outlet_y, outlet_z)
    % Creates a 3D animation of the fluid flow through the grid
    filename = '3D_fluid_flow.gif';
    h = figure('Position', [100 100 800 600], 'Color', 'w');
    
    max_conc = max(c_s(:));
    
    % Create grid coordinates using MESHGRID (MATLAB's expected format)
    [X, Y, Z] = meshgrid(1:gridN_x, 1:gridN_y, 1:gridN_z);
    
    for i = 1:length(t)
        % Get concentration data and permute dimensions to match MESHGRID
        conc3d = squeeze(c_s(i,:,:,:));       % Original dimensions [nx, ny, nz]
        conc3d_permuted = permute(conc3d, [2 1 3]); % Now [ny, nx, nz] to match MESHGRID
        
        % Slice positions to visualize flow path
        xslice = [1, gridN_x];       % Inlet/outlet columns
        yslice = round(gridN_y/2);   % Mid-depth slice
        zslice = [1, gridN_z];       % Bottom/top layers
        
        % Plot concentration slices with corrected data orientation
        slice(X, Y, Z, conc3d_permuted, xslice, yslice, zslice);
        hold on;
        
        % Plot inlet points (corrected indexing for permuted data)
        inlet_y_range = inlet_y(1):inlet_y(2);
        [inX, inY, inZ] = meshgrid(inlet_x, inlet_y_range, inlet_z);
        inC = conc3d_permuted(inlet_y_range, inlet_x, inlet_z);
        scatter3(inX(:), inY(:), inZ(:), 100, inC(:), 'filled', 'MarkerEdgeColor', 'k');
        
        % Plot outlet points (corrected indexing for permuted data)
        outlet_y_range = outlet_y(1):outlet_y(2);
        [outX, outY, outZ] = meshgrid(outlet_x, outlet_y_range, outlet_z);
        outC = conc3d_permuted(outlet_y_range, outlet_x, outlet_z);
        scatter3(outX(:), outY(:), outZ(:), 100, outC(:), 'filled', 'MarkerEdgeColor', 'k');
        
        hold off;
        
        % Styling and labels
        shading interp;
        colormap(jet);
        colorbar;
        caxis([0 max_conc]);
        title(sprintf('3D Fluid Flow at t = %.1f s', t(i)));
        xlabel('X (Length)'); ylabel('Y (Depth)'); zlabel('Z (Height)');
        axis tight;
        view(45, 30); 
        grid on;
        set(gca, 'FontSize', 12);
        
        % Save animation frame
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