function [t_all, c_s, s] = simulate_3d_flow_model_with_pulses(...
    nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid,dx,dz)
    
    % Initialize state variables
    num_cells = nx * ny * nz;
    c_s0 = zeros(nx, ny, nz);
    c_s0(1, :, :) = concentrations(1); % Initial concentration
    s0 = s0_grid;
    %Q0 = zeros(nx, ny, nz);
    %R0 = zeros(nx, ny, nz);
    y0 = [c_s0(:); s0(:)]; %Q0(:); R0(:)];
    
    % Setup ODE options
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
    
    % --- START: PREALLOCATION LOGIC ---
    
    % 1. Calculate the total number of points for preallocation
    num_segments = length(t_breaks) - 1;
    num_points_per_segment = 1100;
    % Total points = points from segment 1 + points from all other segments (excluding duplicates)
    total_points = num_points_per_segment + (num_segments - 1) * (num_points_per_segment - 1);
    
    % 2. Preallocate results arrays using zeros()
    t_all = zeros(total_points, 1);
    y_all = zeros(total_points, length(y0));
    
    % 3. Initialize an index to keep track of where to insert data
    last_idx = 0;
    
    % --- END: PREALLOCATION LOGIC ---
    
    % Process each time segment
    for seg = 1:num_segments
        t_start = t_breaks(seg);
        t_end = t_breaks(seg+1);
        c0_seg = concentrations(seg);
        
        % Determine time points for segment
        tspan = linspace(t_start, t_end, num_points_per_segment);
        
        % Run simulation for segment
        [t_seg, y_seg] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, ...
            kon_grid, koff_grid, smax_grid, c0_seg, D_coeff, ru_to_m,dx, dz), tspan, y0, options);
        
        % --- START: MODIFIED RESULT HANDLING ---
        
        if seg == 1
            % For the first segment, add all points
            num_to_add = num_points_per_segment;
            current_indices = (last_idx + 1):(last_idx + num_to_add);
            t_all(current_indices) = t_seg;
            y_all(current_indices, :) = y_seg;
            last_idx = last_idx + num_to_add;
        else
            % For subsequent segments, skip the first point to avoid duplicates
            num_to_add = num_points_per_segment - 1;
            current_indices = (last_idx + 1):(last_idx + num_to_add);
            t_all(current_indices) = t_seg(2:end);
            y_all(current_indices, :) = y_seg(2:end, :);
            last_idx = last_idx + num_to_add;
        end

        % --- END: MODIFIED RESULT HANDLING ---
        
        % Update initial condition for next segment
        if seg < num_segments
            y0 = y_seg(end, :)';
            c_s_end = reshape(y0(1:num_cells), [nx, ny, nz]);
            s_end = reshape(y0(num_cells+1:2*num_cells), [nx, ny, nz]);
            c_s_end(1, :, :) = concentrations(seg+1);
            y0 = [c_s_end(:); s_end(:); y0(2*num_cells+1:end)]; 
        end
    end
    
    % Extract variables
    c_s = reshape(y_all(:, 1:num_cells), [length(t_all), nx, ny, nz]);
    s = reshape(y_all(:, num_cells+1:2*num_cells), [length(t_all), nx, ny, nz]);
    %Q = reshape(y_all(:, 2*num_cells+1:3*num_cells), [length(t_all), nx, ny, nz]);
    %R = reshape(y_all(:, 3*num_cells+1:end), [length(t_all), nx, ny, nz]);
    
    % Compute kernel
    %K = exp(-Q) .* R;
end