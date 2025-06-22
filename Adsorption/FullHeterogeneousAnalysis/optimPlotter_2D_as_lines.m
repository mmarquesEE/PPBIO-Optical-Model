function stop = optimPlotter_2D_as_lines(log_p, optimValues, state, p_true_2D, exp_data, model_config)
    persistent fig_handle plot_data; 
    if strcmp(state, 'init')
        plot_data.p_true_2D = p_true_2D; 
        plot_data.exp_data = exp_data;
        plot_data.num_sites = (model_config.ads_x_range(2) - model_config.ads_x_range(1) + 1) * (model_config.ads_y_range(2) - model_config.ads_y_range(1) + 1);
        fig_handle = updatePlots_2D_as_lines(log_p, optimValues, plot_data);
    elseif strcmp(state, 'iter')
        if ~ishandle(fig_handle), stop = true; return; end
        updatePlots_2D_as_lines(log_p, optimValues, plot_data);
    elseif strcmp(state, 'done')
        % --- CHANGE: This plotter no longer saves. It just closes. ---
        if ishandle(fig_handle)
            close(fig_handle);
        end
    end
    stop = false;
end