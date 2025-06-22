function stop = optimPlotter_2D(log_p, optimValues, state, p_true_2D, exp_data, model_config)
    % This function is the main callback manager. It handles the optimization
    % state and calls the plotting function.
    persistent fig_handle plot_data; 

    switch state
        case 'init'
            % Store static data that doesn't change during optimization
            plot_data.p_true_2D = p_true_2D;
            plot_data.exp_data = exp_data;
            plot_data.model_config = model_config;
            
            % Call updatePlots_2D to create the figure for the first time
            fig_handle = updatePlots_2D(log_p, optimValues, plot_data);

        case 'iter'
            % Check if the user closed the figure
            if ~ishandle(fig_handle)
                stop = true;
                return; % Exit if figure is closed
            end
            % Call updatePlots_2D to refresh the data
            updatePlots_2D(log_p, optimValues, plot_data);
        case 'done'
            if ishandle(fig_handle)
                % Save the final figure state
                save_pub_fig(fig_handle, 'Adsorption/FullHeterogeneousAnalysis/Figures/figure_fit_diagnostics_2D_FULL', 18);
                fprintf('Full 2D diagnostic plot saved.\n');
            end
    end
    stop = false;
end