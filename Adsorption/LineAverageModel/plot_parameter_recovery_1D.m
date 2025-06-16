function plot_parameter_recovery_1D(p_true, p_opt, p_init, ads_y_dim)
    
    % --- Extract parameter groups from the vectors ---
    true_kon = p_true(1:ads_y_dim);
    opt_kon = p_opt(1:ads_y_dim);
    init_kon = p_init(1:ads_y_dim);
    
    true_koff = p_true(ads_y_dim+1:2*ads_y_dim);
    opt_koff = p_opt(ads_y_dim+1:2*ads_y_dim);
    init_koff = p_init(ads_y_dim+1:2*ads_y_dim);
    
    true_smax = p_true(2*ads_y_dim+1:end);
    opt_smax = p_opt(2*ads_y_dim+1:end);
    init_smax = p_init(2*ads_y_dim+1:end);
    % --- NEW: Calculate Mean Absolute Relative Error (MARE) for each parameter ---
    err_kon_mare = mean(abs(true_kon - opt_kon) ./ abs(true_kon)) * 100;
    err_koff_mare = mean(abs(true_koff - opt_koff) ./ abs(true_koff)) * 100;
    err_smax_mare = mean(abs(true_smax - opt_smax) ./ abs(true_smax)) * 100;
    % --- Create figure ---
    fig5 = figure('Position', [100, 100, 1200, 900]);
    
    % --- Plot kon recovery ---
    ax1 = subplot(3,2,1);
    plot(ax1, true_kon, 'ro-', 'LineWidth', 2, 'DisplayName', 'True'); hold on;
    plot(ax1, init_kon, 'bx--', 'LineWidth', 1.5, 'DisplayName', 'Initial Guess');
    plot(ax1, opt_kon, 'g*-', 'LineWidth', 1.5, 'DisplayName', 'Recovered');
    hold off; grid on; legend('Location', 'best');
    ylabel('Value');
    title(sprintf('k_{on} Recovery (MARE: %.2f%%)', err_kon_mare));
    
    ax2 = subplot(3,2,2);
    bar(ax2, 100 * abs(opt_kon - true_kon) ./ abs(true_kon));
    title('Relative Error in k_{on}');
    ylabel('Error (%)'); grid on;

    % --- Plot koff recovery ---
    ax3 = subplot(3,2,3);
    plot(ax3, true_koff, 'ro-', 'LineWidth', 2); hold on;
    plot(ax3, init_koff, 'bx--', 'LineWidth', 1.5);
    plot(ax3, opt_koff, 'g*-', 'LineWidth', 1.5);
    hold off; grid on;
    ylabel('Value');
    title(sprintf('k_{off} Recovery (MARE: %.2f%%)', err_koff_mare));
    
    ax4 = subplot(3,2,4);
    bar(ax4, 100 * abs(opt_koff - true_koff) ./ abs(true_koff));
    title('Relative Error in k_{off}');
    ylabel('Error (%)'); grid on;
    
    % --- Plot smax recovery ---
    ax5 = subplot(3,2,5);
    plot(ax5, true_smax, 'ro-', 'LineWidth', 2); hold on;
    plot(ax5, init_smax, 'bx--', 'LineWidth', 1.5);
    plot(ax5, opt_smax, 'g*-', 'LineWidth', 1.5);
    hold off; grid on;
    xlabel('Line Index'); ylabel('Value');
    title(sprintf('s_{max} Recovery (MARE: %.2f%%)', err_smax_mare));
    
    ax6 = subplot(3,2,6);
    bar(ax6, 100 * abs(opt_smax - true_smax) ./ abs(true_smax));
    title('Relative Error in s_{max}');
    xlabel('Parameter Index'); ylabel('Error (%)'); grid on;
    
%     sgtitle('Parameter Recovery Results with Quantitative Error', 'FontSize', 16, 'FontWeight', 'bold');
    
    % --- SAVE ENTIRE FIGURE ---
    print(fig5, 'Adsorption/LineAverageModel/Figures/figure_parameter_recovery.png', '-dpng', '-r300');
    print(fig5, 'Adsorption/LineAverageModel/Figures/EPS/figure_parameter_recovery.eps', '-depsc');
    fprintf('Full parameter recovery plot saved as figure_parameter_recovery.png and .eps\n');

    % --- SAVE EACH SUBPLOT AS A SEPARATE EPS FILE ---
    fprintf('Saving individual parameter recovery subplots as EPS files...\n');
    base_path = 'Adsorption/LineAverageModel/Figures/EPS/';
    
    save_subplot_as_eps(ax1, [base_path, 'param_recovery_kon_values.eps']);
    save_subplot_as_eps(ax2, [base_path, 'param_recovery_kon_error.eps']);
    save_subplot_as_eps(ax3, [base_path, 'param_recovery_koff_values.eps']);
    save_subplot_as_eps(ax4, [base_path, 'param_recovery_koff_error.eps']);
    save_subplot_as_eps(ax5, [base_path, 'param_recovery_smax_values.eps']);
    save_subplot_as_eps(ax6, [base_path, 'param_recovery_smax_error.eps']);
    
    fprintf('Finished saving all individual recovery subplots.\n');
end