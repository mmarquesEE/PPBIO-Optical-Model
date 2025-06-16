function stop = optimPlotter_1D(log_params, optimValues, state, true_params_1D, exp_data, ads_y_dim)
    persistent handles; 
    switch state
        case 'init'
            handles.fig = figure('Name', '1D Optimization Progress & Diagnostics', 'Position', [50, 50, 1800, 900]);
            handles.ax_kon = subplot(3, 4, 1); handles.ax_koff = subplot(3, 4, 2); handles.ax_smax = subplot(3, 4, 3);
            handles.ax_fit1 = subplot(3, 4, 5); handles.ax_fit2 = subplot(3, 4, 6); handles.ax_fit3 = subplot(3, 4, 7);
            handles.ax_res1 = subplot(3, 4, 9); handles.ax_res2 = subplot(3, 4, 10); handles.ax_res3 = subplot(3, 4, 11);
            handles.ax_qq = subplot(3, 4, [4, 8, 12]);
            handles.true_params_1D = true_params_1D; handles.exp_data = exp_data; handles.ads_y_dim = ads_y_dim;
            linkaxes([handles.ax_fit1, handles.ax_res1], 'x'); linkaxes([handles.ax_fit2, handles.ax_res2], 'x'); linkaxes([handles.ax_fit3, handles.ax_res3], 'x');
            updatePlots(log_params, optimValues, handles);
        case 'iter'
            if ishandle(handles.fig), updatePlots(log_params, optimValues, handles); else, stop = true; end
        case 'done'
            if ishandle(handles.fig)
%                 sgtitle(handles.fig, 'Optimization Finished!', 'FontSize', 16, 'FontWeight', 'bold');
                % --- SAVE FIGURE ---
                print(handles.fig, 'Adsorption/LineAverageModel/Figures/figure_fit_diagnostics', '-dpng', '-r300');
                fprintf('Diagnostic plot saved as figure_5_fit_diagnostics.png\n');
            end
    end
    stop = false;
end