function updatePlots(current_log_params, optimVals, h)
    % Updates all diagnostic plots during optimization, with errors in the legend.
    
    % --- Parameter Plotting Section (no changes here) ---
    current_params_linear = 10.^current_log_params; ads_y_dim = h.ads_y_dim;
    true_params = h.true_params_1D;
    true_kon = true_params(1:ads_y_dim); true_koff = true_params(ads_y_dim+1:2*ads_y_dim); true_smax = true_params(2*ads_y_dim+1:end);
    opt_kon = current_params_linear(1:ads_y_dim); opt_koff = current_params_linear(ads_y_dim+1:2*ads_y_dim); opt_smax = current_params_linear(2*ads_y_dim+1:end);
    cla(h.ax_kon); plot(h.ax_kon, true_kon, 'ro-'); hold(h.ax_kon, 'on'); plot(h.ax_kon, opt_kon, 'g*-'); hold(h.ax_kon, 'off'); title(h.ax_kon, sprintf('k_{on} (Iter: %d)', optimVals.iteration)); legend(h.ax_kon,{'True','Current'},'Location','best'); grid(h.ax_kon,'on'); xlabel(h.ax_kon, 'Line Index');
    cla(h.ax_koff); plot(h.ax_koff, true_koff, 'ro-'); hold(h.ax_koff, 'on'); plot(h.ax_koff, opt_koff, 'g*-'); hold(h.ax_koff, 'off'); title(h.ax_koff, sprintf('k_{off} (F-count: %d)', optimVals.funccount)); grid(h.ax_koff,'on'); xlabel(h.ax_koff, 'Line Index');
    cla(h.ax_smax); plot(h.ax_smax, true_smax, 'ro-'); hold(h.ax_smax, 'on'); plot(h.ax_smax, opt_smax, 'g*-'); hold(h.ax_smax, 'off'); title(h.ax_smax, sprintf('s_{max} (Residual: %.2e)', optimVals.resnorm)); grid(h.ax_smax,'on'); xlabel(h.ax_smax, 'Line Index');
    
    % --- Sensorgram and Residual Plotting Section ---
    data_struct = h.exp_data{1}; t_exp = data_struct.time; exp_data_matrix = data_struct.signals;
    
    % The 'residual' from the optimizer is defined as (model - data)
    residual_total = optimVals.residual(1:numel(exp_data_matrix));
    residual_matrix = reshape(residual_total, size(exp_data_matrix));
    
    % --- THE FIX IS HERE: Change minus to plus ---
    % To get the smooth model fit, we do: model = data + (model - data)
    s_sim_matrix = exp_data_matrix + residual_matrix;
    
    % --- The rest of the function continues as before ---
    lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);
    ax_fits = [h.ax_fit1, h.ax_fit2, h.ax_fit3]; ax_ress = [h.ax_res1, h.ax_res2, h.ax_res3];
    
    for i = 1:length(lines_to_plot)
        line_idx = lines_to_plot(i); ax_fit = ax_fits(i); ax_res = ax_ress(i);
        
        % Calculate error metrics for the current line's fit
        [~, error_abs] = calculate_validation_error(exp_data_matrix(:, line_idx), s_sim_matrix(:, line_idx));
        fit_legend_text = sprintf('Fit, E_A=%.1e', error_abs);
        
        % Plot the Fit
        cla(ax_fit);
        plot(ax_fit, t_exp, exp_data_matrix(:, line_idx), 'b-', 'DisplayName', 'Data');
        hold(ax_fit, 'on');
        plot(ax_fit, t_exp, s_sim_matrix(:, line_idx), 'r--', 'DisplayName', fit_legend_text);
        hold(ax_fit, 'off');
        title(ax_fit, sprintf('Sensorgram Fit for Line %d', line_idx));
        ylabel(ax_fit, 'Response (RU)');
        legend(ax_fit, 'Location', 'best'); grid(ax_fit, 'on'); set(ax_fit, 'XTickLabel', []);
        
        % Plot the Residuals
        cla(ax_res);
        plot(ax_res, t_exp, residual_matrix(:, line_idx), 'k.', 'MarkerSize', 4);
        hold(ax_res, 'on');
        yline(ax_res, 0, 'r--', 'LineWidth', 1);
        hold(ax_res, 'off');
        title(ax_res, 'Residual (Data - Fit)'); xlabel(ax_res, 'Time (s)'); ylabel(ax_res, 'Error'); grid(ax_res, 'on'); xlim(ax_res, [0, t_exp(end)]);
    end
    
    % --- Q-Q Plot Section ---
    cla(h.ax_qq); qqplot(h.ax_qq, residual_total);
    title(h.ax_qq, 'Q-Q Plot of Residuals'); xlabel(h.ax_qq, 'Standard Normal Quantiles'); ylabel(h.ax_qq, 'Residual Quantiles'); grid(h.ax_qq, 'on');
    
    drawnow;
end