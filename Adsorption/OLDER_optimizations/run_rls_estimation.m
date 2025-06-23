function [theta, P, theta_history] = run_rls_estimation(data, file_name_str, theta_initial, P_initial, lambda_, step_time)
    % --- Prepare data ---
    y = data.RefractiveIndex/1000;
    time = data.Time / 1000;
    dt = mean(diff(time));
    
    % --- Parse concentration from filename ---
    C_str = regexp(file_name_str, '\d+[._]\d*', 'match');
    if isempty(C_str)
        error('Could not parse concentration from filename: %s.', file_name_str);
    end
    concentration_str = strrep(C_str{1}, '_', '.');
    concentration_val = str2double(concentration_str);

    % --- Create the time-varying input vector u(t) ---
    u = zeros(size(time));
    u(time >= step_time) = concentration_val;

    % --- Initialize RLS ---
    theta = theta_initial;
    P = P_initial;
    n_points = length(y);
    theta_history = zeros(n_points, 3);
    
    % --- Main RLS Loop ---
    for k = 1:(n_points - 1)
        phi = [y(k); u(k); u(k)*y(k)];
        y_obs = y(k+1);
        
        K = (P * phi) / (lambda_ + phi' * P * phi);
        theta = theta + K * (y_obs - phi' * theta);
        
        % --- NEW: Apply constraints to keep parameters physical ---
        % theta(1) must be < 1 for kd to be positive
        % theta(2) must be > 0 (assuming ka and N_max are positive)
        % theta(3) must be < 0 for ka to be positive
        theta(1) = min(theta(1), 1 - 1e-9);
        theta(2) = max(theta(2), 1e-9);
        theta(3) = min(theta(3), -1e-9);
        
        P = (P - K * phi' * P) / lambda_;
        theta_history(k, :) = theta';
    end
    theta_history(end, :) = theta';
end