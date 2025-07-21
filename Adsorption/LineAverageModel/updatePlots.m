function fig = updatePlots(current_log_params, optimVals, plot_data)
    
    % --- Find the figure by its Tag. If it doesn't exist, create it. ---
    figTag = 'OptimDashboardFigure';
    fig = findobj('Type', 'figure', 'Tag', figTag);
    background_color = [217/255, 217/255, 217/255];
    if isempty(fig)
        % This block runs only ONCE, during the 'init' phase
        fig = figure('Name', '1D Optimization Progress & Diagnostics', ...
                     'Units', 'centimeters', 'Position', [5, 5, 25, 15], ...
                     'Tag', figTag, 'NumberTitle', 'off', 'Color', background_color);
        tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
        % Just create the tiles. We will find their handles later.
        nexttile(1); nexttile(2); nexttile(3);
        nexttile(5); nexttile(6); nexttile(7);
        nexttile(9); nexttile(10); nexttile(11);
        nexttile(4, [3, 1]); % The spanned Q-Q plot
    end
    
    % --- GET ALL AXES HANDLES (This runs EVERY time) ---
    all_axes = findobj(fig.Children, 'Type', 'Axes');
    % Handles are returned in Last-In, First-Out order (reversed creation order)
    ax_qq   = all_axes(1);
    ax_res3 = all_axes(2); ax_res2 = all_axes(3); ax_res1 = all_axes(4);
    ax_fit3 = all_axes(5); ax_fit2 = all_axes(6); ax_fit1 = all_axes(7);
    ax_smax = all_axes(8); ax_koff = all_axes(9); ax_kon  = all_axes(10);
    
    % --- Define Styles ---
    font_size = 8;
    
    % --- Get static data from the plot_data struct ---
    true_params = plot_data.true_params_1D;
    ads_y_dim = plot_data.ads_y_dim;
    
    % --- Parameter Plotting Section ---
    current_params_linear = 10.^current_log_params;
    true_kon = true_params(1:ads_y_dim); true_koff = true_params(ads_y_dim+1:2*ads_y_dim); true_smax = true_params(2*ads_y_dim+1:end);
    opt_kon = current_params_linear(1:ads_y_dim); opt_koff = current_params_linear(ads_y_dim+1:2*ads_y_dim); opt_smax = current_params_linear(2*ads_y_dim+1:end);
    
    cla(ax_kon); plot(ax_kon, true_kon, 'ro-'); hold(ax_kon, 'on'); plot(ax_kon, opt_kon, 'g*-'); hold(ax_kon, 'off'); title(ax_kon, sprintf('k_{on} (Iter: %d)', optimVals.iteration)); legend(ax_kon,{'True','Current'},'Location','best'); grid(ax_kon,'on'); xlabel(ax_kon, 'Line Index'); set(ax_kon, 'FontSize', font_size-1);
    cla(ax_koff); plot(ax_koff, true_koff, 'ro-'); hold(ax_koff, 'on'); plot(ax_koff, opt_koff, 'g*-'); hold(ax_koff, 'off'); title(ax_koff, sprintf('k_{off} (F-count: %d)', optimVals.funccount)); grid(ax_koff,'on'); xlabel(ax_koff, 'Line Index'); set(ax_koff, 'FontSize', font_size-1);
    cla(ax_smax); plot(ax_smax, true_smax, 'ro-'); hold(ax_smax, 'on'); plot(ax_smax, opt_smax, 'g*-'); hold(ax_smax, 'off'); title(ax_smax, sprintf('s_{max} (Res: %.1e)', optimVals.resnorm)); grid(ax_smax,'on'); xlabel(ax_smax, 'Line Index'); set(ax_smax, 'FontSize', font_size-1);
    
    % The 'residual' field may not exist on the very first call. Check for it.
    if isfield(optimVals, 'residual') && ~isempty(optimVals.residual)
        % --- Sensorgram and Residual Plotting Section ---
        data_struct = plot_data.exp_data{1}; t_exp = data_struct.time; exp_data_matrix = data_struct.signals;
        residual_total = optimVals.residual(1:numel(exp_data_matrix));
        residual_matrix = reshape(residual_total, size(exp_data_matrix));
        s_sim_matrix = exp_data_matrix - residual_matrix; % Corrected: sim = data - residual
        lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);
        ax_fits = [ax_fit1, ax_fit2, ax_fit3]; 
        ax_ress = [ax_res1, ax_res2, ax_res3];
        
        for i = 1:length(lines_to_plot)
            line_idx = lines_to_plot(i); ax_fit = ax_fits(i); ax_res = ax_ress(i);
            
            cla(ax_fit); plot(ax_fit, t_exp, s_sim_matrix(:, line_idx), 'r--'); hold(ax_fit, 'on'); plot(ax_fit, t_exp, exp_data_matrix(:, line_idx), 'b-'); hold(ax_fit, 'off');
            title(ax_fit, sprintf('Fit (Line %d)', line_idx)); ylabel(ax_fit, 'Response (RU)'); grid(ax_fit, 'on'); set(ax_fit, 'XTickLabel', []); set(ax_fit, 'FontSize', font_size-1);
            legend(ax_fit, {'Data', 'Fit'}, 'Location', 'best')

            cla(ax_res); plot(ax_res, t_exp, residual_matrix(:, line_idx), 'k.', 'MarkerSize', 4); hold(ax_res, 'on'); yline(ax_res, 0, 'r--'); hold(ax_res, 'off');
            title(ax_res, 'Residuals'); xlabel(ax_res, 'Time (s)'); ylabel(ax_res, 'Error'); grid(ax_res, 'on'); xlim(ax_res, [0, t_exp(end)]); set(ax_res, 'FontSize', font_size-1);
        end
        
        % --- Q-Q Plot Section ---
        cla(ax_qq); qqplot(ax_qq, residual_total);
        title(ax_qq, 'Q-Q Plot of Residuals'); xlabel(ax_qq, 'Std. Normal Quantiles'); ylabel(ax_qq, 'Residual Quantiles'); grid(ax_qq, 'on'); set(ax_qq, 'FontSize', font_size-1);
    end
    
    drawnow; % Force the figure to update in the event loop
end