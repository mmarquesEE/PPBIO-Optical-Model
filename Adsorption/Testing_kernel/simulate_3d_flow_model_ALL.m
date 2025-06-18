function [t, c_s, s, K, Q] = simulate_3d_flow_model_ALL(nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, c0_assoc, c0_diss, t_assoc, t_total, D_coeff, ru_to_m, s0_grid,dx,dz)
    % Initialize concentrations and auxiliary variables
    c_s = zeros(nx, ny, nz);
    c_s(1, :, :) = c0_assoc;
    s = s0_grid;  % Use provided initial condition
    Q = zeros(nx, ny, nz);
    R = zeros(nx, ny, nz);
    y0 = [c_s(:); s(:); Q(:); R(:)];

    % Time parameters
    tspan_assoc = linspace(0, t_assoc, 1000);
    tspan_diss = linspace(t_assoc, t_total, 1000);
    
    % Solve ODE for association phase
    options = odeset('RelTol',1e-5,'AbsTol',1e-7);
    [t_assoc, y_assoc] = ode15s(@(t,y) ode_system_full(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0_assoc, D_coeff, ru_to_m,dx,dz), tspan_assoc, y0, options);
    
    % Reset for dissociation phase
    y_end_assoc = y_assoc(end,:)';
    num_cells = nx*ny*nz;
    c_s_end = reshape(y_end_assoc(1:num_cells), [nx, ny, nz]);
    s_end = reshape(y_end_assoc(num_cells+1:2*num_cells), [nx, ny, nz]);
    c_s_end(1, :, :) = c0_diss;
    y0_diss = [c_s_end(:); s_end(:); y_end_assoc(2*num_cells+1:end)]; % Preserve Q/R
    
    % Solve ODE for dissociation phase
    [t_diss, y_diss] = ode15s(@(t,y) ode_system_full(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0_diss, D_coeff, ru_to_m,dx,dz), tspan_diss, y0_diss, options);
    
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