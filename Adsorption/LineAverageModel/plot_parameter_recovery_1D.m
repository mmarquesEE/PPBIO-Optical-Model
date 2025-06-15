function plot_parameter_recovery_1D(p_true, p_opt, p_init, ads_y_dim)
    
    % Extract parameter groups
    true_kon = p_true(1:ads_y_dim);
    true_koff = p_true(ads_y_dim+1:2*ads_y_dim);
    true_smax = p_true(2*ads_y_dim+1:end);
    
    opt_kon = p_opt(1:ads_y_dim);
    opt_koff = p_opt(ads_y_dim+1:2*ads_y_dim);
    opt_smax = p_opt(2*ads_y_dim+1:end);

    init_kon = p_init(1:ads_y_dim);
    init_koff = p_init(ads_y_dim+1:2*ads_y_dim);
    init_smax = p_init(2*ads_y_dim+1:end);
    
    figure('Position', [100, 100, 800, 900]);
    
    % Plot kon recovery
    subplot(3,1,1);
    plot(true_kon, 'ro-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_kon, 'bx--', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_kon, 'g*-', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Line-by-Line Binding Rate (k_{on}) Recovery');
    legend('True (Line Avg)', 'Initial Guess', 'Recovered');
    ylabel('Value'); xlabel('Line Index (j)');
    grid on;

    % Plot koff recovery
    subplot(3,1,2);
    plot(true_koff, 'ro-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_koff, 'bx--', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_koff, 'g*-', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Line-by-Line Unbinding Rate (k_{off}) Recovery');
    ylabel('Value'); xlabel('Line Index (j)');
    grid on;
    
    % Plot smax recovery
    subplot(3,1,3);
    plot(true_smax, 'ro-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_smax, 'bx--', 'MarkerSize', 8, 'LineWidth' ,1.5);
    plot(opt_smax, 'g*-', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Line-by-Line Maximum Binding (s_{max}) Recovery');
    ylabel('Value'); xlabel('Line Index (j)');
    grid on;
    
    sgtitle('Parameter Recovery Results for 1D Heterogeneous Model');
end