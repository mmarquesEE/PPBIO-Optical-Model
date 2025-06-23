function plot_model_validation(data_dict, file_names, ka, kd, N_max, step_time)
    % Simulates and plots the model response against actual data.
    
    n_experiments = length(file_names);
    figure('Name', 'Model Validation', 'Position', [150 150 900 600], 'Color', 'w');
    colors = lines(n_experiments);
    legend_entries = {};
    
    for i = 1:n_experiments
        file_name_str = file_names{i};
        data = data_dict.(file_name_str);
        
        % --- Prepare data for simulation ---
        y_actual = data.RefractiveIndex / 1000;
        time = data.Time / 1000;
        dt = mean(diff(time));
        
        % --- CRUCIAL FIX: Parse concentration from the passed filename string ---
        C_str = regexp(file_name_str, '\d+[._]\d*', 'match');
        if isempty(C_str)
            warning('Could not parse concentration from filename %s for plotting. Skipping.', file_name_str);
            continue;
        end
        concentration_str = strrep(C_str{1}, '_', '.');
        concentration_val = str2double(concentration_str);
        u = zeros(size(time));
        u(time >= step_time) = concentration_val;
        
        % --- Simulate the model response ---
        y_hat = zeros(size(y_actual));
        y_hat(1) = y_actual(1);
        
        for k = 1:(length(y_actual) - 1)
            dydt = ka * u(k) * (N_max - y_hat(k)) - kd * y_hat(k);
            y_hat(k+1) = y_hat(k) + dydt * dt;
        end
        
        % --- Plotting ---
        plot(time, y_actual * 1e3, 'Color', colors(i,:), 'LineWidth', 2);
        hold on;
        plot(time, y_hat * 1e3, '--', 'Color', colors(i,:), 'LineWidth', 2);
        
        legend_entries{end+1} = sprintf('Actual %.2fM', concentration_val);
        legend_entries{end+1} = sprintf('Estimated %.2fM', concentration_val);
    end
    
    title('Model Validation: Actual Data vs. Estimated Response');
    xlabel('Time (s)');
    ylabel('\Delta n_s^{eff} [mRIU]');
    legend(legend_entries, 'Location', 'best');
    grid on;
    hold off;
end