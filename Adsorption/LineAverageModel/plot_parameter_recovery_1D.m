function plot_parameter_recovery_1D(p_true, p_opt, p_init, ads_y_dim, y_errors_neg, y_errors_pos)
    % --- Define Publication Style Parameters ---
    base_font_size = 8;
    line_width = 1.2;
    marker_size = 12;
    error_cap_size = 4;

    % --- Extract parameter groups from the vectors (code is unchanged) ---
    true_kon = p_true(1:ads_y_dim); opt_kon = p_opt(1:ads_y_dim); init_kon = p_init(1:ads_y_dim);
    true_koff = p_true(ads_y_dim+1:2*ads_y_dim); opt_koff = p_opt(ads_y_dim+1:2*ads_y_dim); init_koff = p_init(ads_y_dim+1:2*ads_y_dim);
    true_smax = p_true(2*ads_y_dim+1:end); opt_smax = p_opt(2*ads_y_dim+1:end); init_smax = p_init(2*ads_y_dim+1:end);
    y_err_neg_kon = y_errors_neg(1:ads_y_dim); y_err_pos_kon = y_errors_pos(1:ads_y_dim);
    y_err_neg_koff = y_errors_neg(ads_y_dim+1:2*ads_y_dim); y_err_pos_koff = y_errors_pos(ads_y_dim+1:2*ads_y_dim);
    y_err_neg_smax = y_errors_neg(2*ads_y_dim+1:end); y_err_pos_smax = y_errors_pos(2*ads_y_dim+1:end);
    
    % --- Calculate MARE (code is unchanged) ---
    err_kon_mare = mean(abs(true_kon - opt_kon) ./ abs(true_kon)) * 100;
    err_koff_mare = mean(abs(true_koff - opt_koff) ./ abs(true_koff)) * 100;
    err_smax_mare = mean(abs(true_smax - opt_smax) ./ abs(true_smax)) * 100;
    
    line_indices = 1:ads_y_dim;

    % =========================================================================
    % --- FIGURE 1: The "Simplified" Plot for the Main Paper ---
    % =========================================================================
    fprintf('\nGenerating simplified 3x1 recovery plot for publication (8.4cm wide)...\n');
    fig_pub = figure('Name', 'Simplified Parameter Recovery');
    t_pub = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Clean plot for k_on ---
    ax1_pub = nexttile;
    plot(ax1_pub, line_indices, true_kon, 'ro-', 'LineWidth', line_width, 'DisplayName', 'True Value'); hold(ax1_pub, 'on');
    errorbar(ax1_pub, line_indices, opt_kon, y_err_neg_kon, y_err_pos_kon, 'k.', 'MarkerSize', marker_size, 'CapSize', error_cap_size, 'LineWidth', line_width-0.2, 'DisplayName', 'Estimated (95% CI)', 'LineStyle', 'none');
    hold(ax1_pub, 'off'); grid on; box on; legend('Location', 'best', 'FontSize', base_font_size-1);
    ylabel(ax1_pub, 'Value', 'FontSize', base_font_size); title(ax1_pub, 'k_{on} Recovery', 'FontSize', base_font_size);
    xlim(ax1_pub, [0, ads_y_dim + 1]); set(ax1_pub, 'FontSize', base_font_size-1);

    % --- Clean plot for k_off ---
    ax2_pub = nexttile;
    plot(ax2_pub, line_indices, true_koff, 'ro-', 'LineWidth', line_width, 'DisplayName', 'True Value'); hold(ax2_pub, 'on');
    errorbar(ax2_pub, line_indices, opt_koff, y_err_neg_koff, y_err_pos_koff, 'k.', 'MarkerSize', marker_size, 'CapSize', error_cap_size, 'LineWidth', line_width-0.2, 'DisplayName', 'Estimated (95% CI)', 'LineStyle', 'none');
    hold(ax2_pub, 'off'); grid on; box on; legend('Location', 'best', 'FontSize', base_font_size-1);
    ylabel(ax2_pub, 'Value', 'FontSize', base_font_size); title(ax2_pub, 'k_{off} Recovery', 'FontSize', base_font_size);
    xlim(ax2_pub, [0, ads_y_dim + 1]); set(ax2_pub, 'FontSize', base_font_size-1);
    
    % --- Clean plot for s_max ---
    ax3_pub = nexttile;
    plot(ax3_pub, line_indices, true_smax, 'ro-', 'LineWidth', line_width, 'DisplayName', 'True Value'); hold(ax3_pub, 'on');
    errorbar(ax3_pub, line_indices, opt_smax, y_err_neg_smax, y_err_pos_smax, 'k.', 'MarkerSize', marker_size, 'CapSize', error_cap_size, 'LineWidth', line_width-0.2, 'DisplayName', 'Estimated (95% CI)', 'LineStyle', 'none');
    hold(ax3_pub, 'off'); grid on; box on; legend('Location', 'best', 'FontSize', base_font_size-1);
    xlabel(ax3_pub, 'Line Index', 'FontSize', base_font_size); ylabel(ax3_pub, 'Value', 'FontSize', base_font_size);
    title(ax3_pub, 's_{max} Recovery', 'FontSize', base_font_size);
    xlim(ax3_pub, [0, ads_y_dim + 1]); set(ax3_pub, 'FontSize', base_font_size-1);
    
    title(t_pub, 'Final Parameter Recovery: True vs. Estimated', 'FontSize', base_font_size+1, 'FontWeight', 'bold');
    
    % --- Save the simplified, publication-ready figure ---
    save_pub_fig(fig_pub, 'Adsorption/LineAverageModel/Figures/figure_parameter_recovery_simplified_pub', 8.4);
    close(fig_pub);
    
    % =========================================================================
    % --- FIGURE 2: The "Detailed" Plot for Supplementary Materials ---
    % =========================================================================
    fprintf('\nGenerating detailed 3x2 recovery plot for supplementary materials (18cm wide)...\n');
    fig_supp = figure('Name', 'Detailed Parameter Recovery');
    t_supp = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Plot kon recovery ---
    ax1_s = nexttile;
    plot(ax1_s, line_indices, true_kon, 'ro-', 'LineWidth', line_width); hold(ax1_s, 'on');
    plot(ax1_s, line_indices, init_kon, 'bx--', 'LineWidth', line_width);
    errorbar(ax1_s, line_indices, opt_kon, y_err_neg_kon, y_err_pos_kon, 'g*', 'LineWidth', line_width, 'LineStyle', 'none');
    hold(ax1_s, 'off'); grid on; legend(ax1_s, {'True', 'Initial', 'Recovered (95% CI)'},'Location', 'best', 'FontSize', base_font_size-1);
    ylabel(ax1_s, 'Value', 'FontSize', base_font_size); title(ax1_s, sprintf('k_{on} Recovery (MARE: %.2f%%)', err_kon_mare), 'FontSize', base_font_size);
    set(ax1_s, 'FontSize', base_font_size-1);
    
    ax2_s = nexttile;
    bar(ax2_s, line_indices, 100 * abs(opt_kon - true_kon) ./ abs(true_kon));
    title(ax2_s, 'Relative Error in k_{on}', 'FontSize', base_font_size); ylabel(ax2_s, 'Error (%)', 'FontSize', base_font_size); grid on;
    set(ax2_s, 'FontSize', base_font_size-1);
    
    % --- Plot koff recovery & error---
    ax3_s = nexttile;
    plot(ax3_s, line_indices, true_koff, 'ro-', 'LineWidth', line_width); hold on; plot(ax3_s, line_indices, init_koff, 'bx--', 'LineWidth', line_width); errorbar(ax3_s, line_indices, opt_koff, y_err_neg_koff, y_err_pos_koff, 'g*', 'LineWidth', line_width, 'LineStyle', 'none');
    hold off; grid on; ylabel(ax3_s, 'Value', 'FontSize', base_font_size); title(ax3_s, sprintf('k_{off} Recovery (MARE: %.2f%%)', err_koff_mare), 'FontSize', base_font_size);
    set(ax3_s, 'FontSize', base_font_size-1);

    ax4_s = nexttile;
    bar(ax4_s, line_indices, 100 * abs(opt_koff - true_koff) ./ abs(true_koff));
    title(ax4_s, 'Relative Error in k_{off}', 'FontSize', base_font_size); ylabel(ax4_s, 'Error (%)', 'FontSize', base_font_size); grid on;
    set(ax4_s, 'FontSize', base_font_size-1);

    % --- Plot smax recovery & error ---
    ax5_s = nexttile;
    plot(ax5_s, line_indices, true_smax, 'ro-', 'LineWidth', line_width); hold on; plot(ax5_s, line_indices, init_smax, 'bx--', 'LineWidth', line_width); errorbar(ax5_s, line_indices, opt_smax, y_err_neg_smax, y_err_pos_smax, 'g*', 'LineWidth', line_width, 'LineStyle', 'none');
    hold off; grid on; xlabel(ax5_s, 'Line Index', 'FontSize', base_font_size); ylabel(ax5_s, 'Value', 'FontSize', base_font_size); title(ax5_s, sprintf('s_{max} Recovery (MARE: %.2f%%)', err_smax_mare), 'FontSize', base_font_size);
    set(ax5_s, 'FontSize', base_font_size-1);

    ax6_s = nexttile;
    bar(ax6_s, line_indices, 100 * abs(opt_smax - true_smax) ./ abs(true_smax));
    title(ax6_s, 'Relative Error in s_{max}', 'FontSize', base_font_size); xlabel(ax6_s, 'Line Index', 'FontSize', base_font_size); ylabel(ax6_s, 'Error (%)', 'FontSize', base_font_size); grid on;
    set(ax6_s, 'FontSize', base_font_size-1);
    
    % --- Save the detailed, supplementary figure ---
    save_pub_fig(fig_supp, 'Adsorption/LineAverageModel/Figures/figure_parameter_recovery_detailed_pub', 18);
    close(fig_supp);

    % The old code for saving individual subplots is no longer needed.
    fprintf('\nFinished generating parameter recovery plots.\n');
end