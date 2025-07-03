function stop = optimPlotter_1D(log_params, optimValues, state, true_params_1D, exp_data, ads_y_dim)
    persistent fig_handle plot_data; % Only persist the main figure handle and static data

    switch state
        case 'init'
            % Create a struct to hold static data that doesn't change
            plot_data.true_params_1D = true_params_1D;
            plot_data.exp_data = exp_data;
            plot_data.ads_y_dim = ads_y_dim;
            
            % Call updatePlots to create and draw the figure for the first time
            fig_handle = updatePlots(log_params, optimValues, plot_data);

        case 'iter'
            % Before each update, check if the user closed the figure
            if ~ishandle(fig_handle)
                stop = true; % Stop the optimizer
                return;
            end
            % Call updatePlots to refresh the data on the existing figure
            updatePlots(log_params, optimValues, plot_data);

        case 'done'
            if ishandle(fig_handle)
                % On the final call, save the figure
                save_pub_fig(fig_handle, 'Adsorption/LineAverageModel/Figures/figure_fit_diagnostics_FULL', 18);
                fprintf('Full diagnostic plot saved for supplementary materials.\n');
            end
    end
    stop = false; % Default to not stopping the optimization
end