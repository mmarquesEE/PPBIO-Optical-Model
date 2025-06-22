function fig = updatePlots_2D_as_lines(log_p, optimVals, plot_data)
    figTag = 'OptimDashboardFigure_2D_Lines';
    fig = findobj('Type', 'figure', 'Tag', figTag);
    if isempty(fig)
        fig = figure('Name', '2D Optimization Progress (Line Plots)', 'Units', 'centimeters', 'Position', [5, 5, 25, 15], 'Tag', figTag, 'NumberTitle', 'off');
        tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
        nexttile(1); nexttile(2); nexttile(3); nexttile(5); nexttile(6); nexttile(7);
        nexttile(9); nexttile(10); nexttile(11); nexttile(4, [3, 1]); 
    end
    all_axes = findobj(fig.Children, 'Type', 'Axes');
    ax_qq   = all_axes(1); ax_res3 = all_axes(2); ax_res2 = all_axes(3); ax_res1 = all_axes(4);
    ax_fit3 = all_axes(5); ax_fit2 = all_axes(6); ax_fit1 = all_axes(7);
    ax_smax = all_axes(8); ax_koff = all_axes(9); ax_kon  = all_axes(10);
    font_size = 8; line_width = 1.2; marker_size = 4;
    
    true_params = plot_data.p_true_2D; num_sites = plot_data.num_sites;
    site_indices = 1:num_sites;

    current_params_linear = 10.^log_p;
    true_kon = true_params(1:num_sites); opt_kon = current_params_linear(1:num_sites);
    true_koff = true_params(num_sites+1:2*num_sites); opt_koff = current_params_linear(num_sites+1:2*num_sites);
    true_smax = true_params(2*num_sites+1:end); opt_smax = current_params_linear(2*num_sites+1:end);
    
    cla(ax_kon); plot(ax_kon, site_indices, true_kon, 'ro-'); hold(ax_kon, 'on'); plot(ax_kon, site_indices, opt_kon, 'g.-'); hold(ax_kon, 'off'); title(ax_kon, 'k_{on}'); grid(ax_kon,'on'); xlabel(ax_kon, 'Site Index');
    cla(ax_koff); plot(ax_koff, site_indices, true_koff, 'ro-'); hold(ax_koff, 'on'); plot(ax_koff, site_indices, opt_koff, 'g.-'); hold(ax_koff, 'off'); title(ax_koff, 'k_{off}'); grid(ax_koff,'on'); xlabel(ax_koff, 'Site Index');
    cla(ax_smax); plot(ax_smax, site_indices, true_smax, 'ro-'); hold(ax_smax, 'on'); plot(ax_smax, site_indices, opt_smax, 'g.-'); hold(ax_smax, 'off'); title(ax_smax, 's_{max}'); grid(ax_smax,'on'); xlabel(ax_smax, 'Site Index');

    if isfield(optimVals, 'residual') && ~isempty(optimVals.residual)
        data_struct = plot_data.exp_data{1}; t_exp = data_struct.time; exp_data_matrix = data_struct.signals;
        residual_total = optimVals.residual(1:numel(exp_data_matrix));
        residual_matrix = reshape(residual_total, size(exp_data_matrix));
        s_sim_matrix = exp_data_matrix + residual_matrix;
        lines_to_plot = unique([1, round(size(exp_data_matrix,2)/2), size(exp_data_matrix,2)]);
        ax_fits = [ax_fit1, ax_fit2, ax_fit3]; ax_ress = [ax_res1, ax_res2, ax_res3];
        
        for i = 1:length(lines_to_plot)
            line_idx = lines_to_plot(i); ax_fit = ax_fits(i); ax_res = ax_ress(i);
            cla(ax_fit); plot(ax_fit, t_exp, exp_data_matrix(:, line_idx), 'b-'); hold(ax_fit, 'on'); plot(ax_fit, t_exp, s_sim_matrix(:, line_idx), 'r--'); hold(ax_fit, 'off');
            title(ax_fit, sprintf('Fit (Line %d)', line_idx)); ylabel(ax_fit, 'Response (RU)'); grid(ax_fit, 'on'); set(ax_fit, 'XTickLabel', []);
            cla(ax_res); plot(ax_res, t_exp, residual_matrix(:, line_idx), 'k.', 'MarkerSize', marker_size); hold(ax_res, 'on'); yline(ax_res, 0, 'r--'); hold(ax_res, 'off');
            title(ax_res, 'Residuals'); xlabel(ax_res, 'Time (s)'); ylabel(ax_res, 'Error'); grid(ax_res, 'on'); xlim(ax_res, [0, t_exp(end)]);
        end
        cla(ax_qq); qqplot(ax_qq, residual_total); title(ax_qq, 'Q-Q Plot of Residuals'); grid(ax_qq, 'on');
    end
    drawnow;
end