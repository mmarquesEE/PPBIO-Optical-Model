function plot_final_2D_results_full(p_true_2D, opt_log_params, p_init, resnorm, final_residual, J_opt, output, exp_data, model_config)
    % Creates a final 3x3 summary with parameter CIs, fits, residuals, Q-Q plot,
    % and now includes the initial parameter guess and total iterations.

    fprintf('\nGenerating comprehensive 3x3 results summary figure with confidence intervals...\n');

    % --- Define Publication Style Parameters ---
    target_fig_width_cm = 18;  
    base_font_size = 8;        
    line_width = 1.0;
    marker_size = 10;
    error_cap_size = 0;

    % --- Get Dimensions ---
    num_sites = (model_config.ads_x_range(2) - model_config.ads_x_range(1) + 1) * ...
                (model_config.ads_y_range(2) - model_config.ads_y_range(1) + 1);
    N_params_2D = 3 * num_sites;
    
    % --- Step 1: Calculate Confidence Intervals ---
    dof = numel(final_residual) - N_params_2D; 
    noise_variance_est = resnorm / dof; 
    covariance_matrix_log = noise_variance_est * pinv(full(J_opt' * J_opt));
    param_variances_log = diag(covariance_matrix_log);
    param_stderr_log = sqrt(param_variances_log);
    ci_95_log = 1.96 * param_stderr_log;
    upper_bound_log = opt_log_params + ci_95_log;
    lower_bound_log = opt_log_params - ci_95_log;
    opt_params_2D = 10.^opt_log_params;
    upper_bound_lin = 10.^upper_bound_log;
    lower_bound_lin = 10.^lower_bound_log;
    y_errors_pos = upper_bound_lin - opt_params_2D;
    y_errors_neg = opt_params_2D - lower_bound_lin;
    
    % --- Create the master figure and 3x3 layout ---
    fig = figure('Name', 'Comprehensive Final Results');
    t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- De-vectorize all necessary data ---
    site_indices = 1:num_sites;
    true_kon = p_true_2D(1:num_sites); opt_kon = opt_params_2D(1:num_sites); init_kon = p_init(1:num_sites);
    err_neg_kon = y_errors_neg(1:num_sites); err_pos_kon = y_errors_pos(1:num_sites);
    
    true_koff = p_true_2D(num_sites+1:2*num_sites); opt_koff = opt_params_2D(num_sites+1:2*num_sites); init_koff = p_init(num_sites+1:2*num_sites);
    err_neg_koff = y_errors_neg(num_sites+1:2*num_sites); err_pos_koff = y_errors_pos(num_sites+1:2*num_sites);

    true_smax = p_true_2D(2*num_sites+1:end); opt_smax = opt_params_2D(2*num_sites+1:end); init_smax = p_init(2*num_sites+1:end);
    err_neg_smax = y_errors_neg(2*num_sites+1:end); err_pos_smax = y_errors_pos(2*num_sites+1:end);

    exp_data_matrix = exp_data{1}.signals;
    t_exp = exp_data{1}.time;
    residual_matrix = reshape(final_residual, size(exp_data_matrix));
    s_sim_matrix = exp_data_matrix + residual_matrix;
    lines_to_plot = unique([1, round(size(exp_data_matrix,2)/2), size(exp_data_matrix,2)]);

    % =====================================================================
    % --- ROW 1: Final Parameter Recovery with Confidence Intervals & Initial Guess ---
    % =====================================================================
    ax1 = nexttile; 
    plot(ax1, site_indices, true_kon, 'r-'); hold(ax1, 'on');
    plot(ax1, site_indices, init_kon, 'bx--');
    errorbar(ax1, site_indices, opt_kon, err_neg_kon, err_pos_kon, 'k.', 'MarkerSize', marker_size, 'CapSize', error_cap_size); 
    hold(ax1, 'off'); grid on; box on; title('k_{on} Recovery'); set(ax1, 'XTickLabel', []); ylabel('Value'); 
    legend(ax1, {'True', 'Initial', 'Estimated (95% CI)'}, 'Location','best');

    ax2 = nexttile; 
    plot(ax2, site_indices, true_koff, 'r-'); hold(ax2, 'on');
    plot(ax2, site_indices, init_koff, 'bx--');
    errorbar(ax2, site_indices, opt_koff, err_neg_koff, err_pos_koff, 'k.', 'MarkerSize', marker_size, 'CapSize', error_cap_size); 
    hold(ax2, 'off'); grid on; box on; title('k_{off} Recovery'); set(ax2, 'XTickLabel', []);

    ax3 = nexttile; 
    plot(ax3, site_indices, true_smax, 'r-'); hold(ax3, 'on');
    plot(ax3, site_indices, init_smax, 'bx--');
    errorbar(ax3, site_indices, opt_smax, err_neg_smax, err_pos_smax, 'k.', 'MarkerSize', marker_size, 'CapSize', error_cap_size); 
    hold(ax3, 'off'); grid on; box on; title('s_{max} Recovery'); set(ax3, 'XTickLabel', []);

    % =====================================================================
    % --- ROW 2: Final Sensorgram Fits --- (No changes here)
    % =====================================================================
    for i = 1:length(lines_to_plot)
        ax = nexttile; line_idx = lines_to_plot(i);
        plot(ax, t_exp, exp_data_matrix(:, line_idx), 'b-'); hold(ax, 'on'); 
        plot(ax, t_exp, s_sim_matrix(:, line_idx), 'r--'); hold(ax, 'off');
        grid on; box on; title(sprintf('Fit (Line %d)', line_idx)); set(ax, 'XTickLabel', []);
        if i == 1, ylabel('Response (RU)'); end
    end
    
    % =====================================================================
    % --- ROW 3: Final Residuals and Q-Q Plot --- (No changes here)
    % =====================================================================
    for i = 1:length(lines_to_plot) - 1
        ax = nexttile; line_idx = lines_to_plot(i);
        plot(ax, t_exp, residual_matrix(:, line_idx), 'k.', 'MarkerSize', marker_size-2);
        hold(ax, 'on'); yline(ax, 0, 'r--'); hold(ax, 'off');
        grid on; box on; title('Residuals'); axis(ax,'tight');
        if i == 1, ylabel('Error (RU)'); end
    end
    ax_qq = nexttile;
    qqplot(ax_qq, final_residual);
    grid(ax_qq, 'on'); box(ax_qq, 'on');
    title(ax_qq, 'Q-Q Plot of All Residuals');

    % --- Apply consistent font sizes and add shared labels ---
    set(findobj(fig, 'Type', 'Axes'), 'FontSize', base_font_size-1);
    xlabel(t, 'Time (s) / Site Index', 'FontSize', base_font_size);
    
    % --- FIX: Add iteration count to the main title ---
    num_iter = output.iterations;
%     main_title_str = sprintf('Comprehensive Post-Optimization Results (Iterations: %d)', num_iter);
%     title(t, main_title_str, 'FontSize', base_font_size+2, 'FontWeight', 'bold');
    
    % --- Save the Final Figure ---
    save_pub_fig(fig, 'FullHeterogeneousAnalysis/Figures/figure_final_summary_full_reg', target_fig_width_cm);
    close(fig);
end