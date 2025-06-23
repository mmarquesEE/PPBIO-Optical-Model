function plot_parameter_convergence(theta_history, dt_mean)
    % Plots the evolution of the physical parameters ka, kd, N_max over time.
    ka_hist = -theta_history(:,3) / dt_mean;
    kd_hist = (1 - theta_history(:,1)) / dt_mean;
    N_max_hist = -theta_history(:,2) ./ theta_history(:,3);
    
    figure('Name', 'RLS Parameter Convergence', 'Position', [100 100 800 600], 'Color', 'w');
    sgtitle('Convergence of Physical Parameters');
    
    subplot(3, 1, 1);
    plot(ka_hist, 'b', 'LineWidth', 1.5);
    ylabel('k_a'); grid on;
    
    subplot(3, 1, 2);
    plot(kd_hist, 'g', 'LineWidth', 1.5);
    ylabel('k_d'); grid on;
    
    subplot(3, 1, 3);
    plot(N_max_hist, 'r', 'LineWidth', 1.5);
    ylabel('N_{max}'); xlabel('Iteration (Total Data Points)'); grid on;
end
