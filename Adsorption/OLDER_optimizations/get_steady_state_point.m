function [R, C] = get_steady_state_point(data, file_name_str, step_time)
    % Calculates a representative steady-state point from data
    time = data.Time / 1000;
    y = data.RefractiveIndex / 1000;
    
    % CRUCIAL FIX: Parse concentration from the passed filename string
    C_str = regexp(file_name_str, '\d+[._]\d*', 'match');
    if isempty(C_str)
        error('Could not parse concentration from filename for initialization: %s.', file_name_str);
    end
    concentration_str = strrep(C_str{1}, '_', '.');
    C = str2double(concentration_str);

    % Get the mean of the last few points as the steady state value
    steady_state_indices = find(time > (time(end) - 30)); % Use last 30s
    if isempty(steady_state_indices) || length(steady_state_indices) < 2
        steady_state_indices = round(0.9*length(y)):length(y); % Fallback: last 10%
    end
    R = mean(y(steady_state_indices));
end