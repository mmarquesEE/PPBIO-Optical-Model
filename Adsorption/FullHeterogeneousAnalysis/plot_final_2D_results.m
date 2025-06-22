function plot_final_2D_results(p_true_2D, opt_log_params, p_init, final_residual, output, exp_data, model_config)
    % Creates a final 3x3 summary figure from the optimization results,
    % including the initial guess and total iterations.

    fprintf('\nGenerating final 3x3 results summary figure...\n');

    % --- Define Publication Style Parameters ---
    target_fig_width_cm = 18;  % A 3x3 grid needs a wider format
    base_font_size = 8;        
    line_width = 1.2;
    marker_size = 4;

    % --- Get Dimensions & Indices ---
    num_sites = (model_config.ads_x_range(2) - model_config.ads_x_range(1) + 1) * ...
                (model_config.ads_y_range(2) - model_config.ads_y_range(1) + 1);
    site_indices = 1:num_sites;
    
    % --- Create the master figure and 3x3 layout ---
    fig = figure('Name', 'Final Optimization Results Summary (No QQ)');
    t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Process Final Parameters and Residuals ---
    opt_params_2D = 10.^opt_log_params;
    
    true_kon = p_true_2D(1:num_sites); opt_kon = opt_params_2D(1:num_sites); init_kon = p_init(1:num_sites);
    true_koff = p_true_2D(num_sites+1:2*num_sites); opt_koff = opt_params_2D(num_sites+1:2*num_sites); init_koff = p_init(num_sites+1:2*num_sites);
    true_smax = p_true_2D(2*num_sites+1:end); opt_smax = opt_params_2D(2*num_sites+1:end); init_smax = p_init(2*num_sites+1:end);
    
    exp_data_matrix = exp_data{1}.signals;
    t_exp = exp_data{1}.time;
    residual_matrix = reshape(final_residual, size(exp_data_matrix));
    s_sim_matrix = exp_data_matrix + residual_matrix;
    lines_to_plot = unique([1, round(size(exp_data_matrix,2)/2), size(exp_data_matrix,2)]);

    % =====================================================================
    % --- ROW 1: Final Parameter Recovery with Initial Guess ---
    % =====================================================================
    ax1 = nexttile; 
    plot(ax1, site_indices, true_kon, 'ro-'); hold(ax1, 'on'); 
    plot(ax1, site_indices, init_kon, 'bx--');
    plot(ax1, site_indices, opt_kon, 'k.-'); 
    hold(ax1, 'off'); grid on; box on; title('k_{on} Recovery'); set(ax1, 'XTickLabel', []); ylabel('Value'); 
    legend({'True','Initial','Estimated'}, 'Location','best');
    
    ax2 = nexttile; 
    plot(ax2, site_indices, true_koff, 'ro-'); hold(ax2, 'on'); 
    plot(ax2, site_indices, init_koff, 'bx--');
    plot(ax2, site_indices, opt_koff, 'k.-'); 
    hold(ax2, 'off'); grid on; box on; title('k_{off} Recovery'); set(ax2, 'XTickLabel', []);
    
    ax3 = nexttile; 
    plot(ax3, site_indices, true_smax, 'ro-'); hold(ax3, 'on'); 
    plot(ax3, site_indices, init_smax, 'bx--');
    plot(ax3, site_indices, opt_smax, 'k.-'); 
    hold(ax3, 'off'); grid on; box on; title('s_{max} Recovery'); set(ax3, 'XTickLabel', []);

    % =====================================================================
    % --- ROW 2: Final Sensorgram Fits ---
    % =====================================================================
    for i = 1:length(lines_to_plot)
        ax = nexttile;
        line_idx = lines_to_plot(i);
        plot(ax, t_exp, exp_data_matrix(:, line_idx), 'b-', 'DisplayName', 'Data'); hold(ax, 'on'); 
        plot(ax, t_exp, s_sim_matrix(:, line_idx), 'r--', 'DisplayName', 'Fit'); hold(ax, 'off');
        grid on; box on;
        title(sprintf('Fit (Line %d)', line_idx));
        set(ax, 'XTickLabel', []);
        if i == 1, ylabel('Response (RU)'); end
        if i == length(lines_to_plot), legend('Location','best'); end
    end
    
    % =====================================================================
    % --- ROW 3: Final Residuals ---
    % =====================================================================
    for i = 1:length(lines_to_plot)
        ax = nexttile;
        line_idx = lines_to_plot(i);
        plot(ax, t_exp, residual_matrix(:, line_idx), 'k.', 'MarkerSize', marker_size);
        hold(ax, 'on'); yline(ax, 0, 'r--'); hold(ax, 'off');
        grid on; box on;
        title('Residuals');
        axis(ax,'tight');
        if i == 1, ylabel('Error (RU)'); end
    end
    
    % --- Apply consistent font sizes and add shared labels ---
    set(findobj(fig, 'Type', 'Axes'), 'FontSize', base_font_size-1);
    xlabel(t, 'Time (s) / Site Index', 'FontSize', base_font_size);
    
    % --- FIX: Add iteration count to the main title ---
    num_iter = output.iterations;
%     main_title_str = sprintf('Final Post-Optimization Results (Iterations: %d)', num_iter);
%     title(t, main_title_str, 'FontSize', base_font_size+2, 'FontWeight', 'bold');
    
    % --- Save the Final Figure ---
    save_pub_fig(fig, 'FullHeterogeneousAnalysis/Figures/figure_final_summary_no_qq', target_fig_width_cm);
    close(fig);
end