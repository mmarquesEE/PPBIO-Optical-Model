function exp_settings = generate_experiments(M, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2)
    % Generates M experiments with orthogonal concentration profiles and flow velocities
    %
    % Inputs:
    %   M - Number of experiments
    %   base_max_velocity - Reference flow velocity (e.g., 8.3)
    %   T1, T2, T3 - Fixed pulse times
    %   t_total - Total experiment duration
    %   c_diss - Dissociation concentration
    %   c1, c2 - Base concentrations
    
    exp_settings = struct(...
        'pulse_times', {}, ...
        'pulse_concs', {}, ...
        't_total', {}, ...
        'max_velocity', {}, ...
        'c_diss', {} ...
    );
    
    % Generate orthogonal concentration pairs using polar coordinates
    angles = linspace(0, pi/2, M);  % Cover quadrant for positive concentrations
    factors = linspace(0.3, 3, M);  % Concentration scaling factors
    
    for i = 1:M
        if M ==1
            M=2;
        end
        % Create orthogonal concentration profiles
        conc_factor1 = factors(ceil(i/2)) * cos(angles(i));
        conc_factor2 = factors(ceil(i/2)) * sin(angles(i));
        
        % Ensure minimum concentration variation
        min_conc = 0.1 * min(c1, c2);
        conc1 = max(c1 * (0.5 + conc_factor1), min_conc);
        conc3 = max(c2 * (0.5 + conc_factor2), min_conc);
        
        % Create velocity profile (logarithmic spacing)
        vel_min = 0.1 * base_max_velocity;
        vel_max = 5.0 * base_max_velocity;
        velocity = exp(log(vel_min) + (i-1)/(M-1) * (log(vel_max) - log(vel_min)));
        
        % Special patterns for every 3rd experiment
        if mod(i,3) == 0
            exp_settings(i).pulse_concs = [conc1, c2, conc3];  % Middle pulse active
        elseif mod(i,4) == 0
            exp_settings(i).pulse_concs = [c1, 0, conc3];      % First pulse fixed
        else
            exp_settings(i).pulse_concs = [conc1, 0, conc3];   % Standard pattern
        end
        
        % Assign common parameters
        exp_settings(i).pulse_times = [T1, T2, T3];
        exp_settings(i).t_total = t_total;
        exp_settings(i).max_velocity = velocity;
        exp_settings(i).c_diss = c_diss;
    end
end