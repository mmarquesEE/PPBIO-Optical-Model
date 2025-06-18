function plot_and_save_parameter_fields(kon_homog, koff_homog, smax_homog, ...
                                        kon_heterog, koff_heterog, smax_heterog, ...
                                        ads_x_range, ads_y_range, ads_layer)
    
    figure('Position', [300, 300, 500, 400]);
    
    % Use tiledlayout for better control over spacing
    tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Row 1: Plot k_on ---
    % Calculate color limits ONLY from the sensible region
    kon_slice_heterog = kon_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_kon = [min(kon_slice_heterog(:)), max(kon_slice_heterog(:))];
    if diff(c_limits_kon) < 1e-9; c_limits_kon(2) = c_limits_kon(1) + 1; end % Robustness check

    % Homogeneous k_on
    nexttile;
    plot_grid_with_black_background(kon_homog(:,:,ads_layer)', c_limits_kon);
    title('Homogeneous');
    ylabel('k_{on}', 'FontSize', 10);
    set(gca, 'XTickLabel', []); 
    
    % Heterogeneous k_on
    nexttile;
    plot_grid_with_black_background(kon_heterog(:,:,ads_layer)', c_limits_kon);
    title('Heterogeneous');
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    colorbar;

    % --- Row 2: Plot k_off ---
    koff_slice_heterog = koff_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_koff = [min(koff_slice_heterog(:)), max(koff_slice_heterog(:))];
    if diff(c_limits_koff) < 1e-9; c_limits_koff(2) = c_limits_koff(1) + 1; end

    % Homogeneous k_off
    nexttile;
    plot_grid_with_black_background(koff_homog(:,:,ads_layer)', c_limits_koff);
    ylabel('k_{off}', 'FontSize', 10);
    set(gca, 'XTickLabel', []); 

    % Heterogeneous k_off
    nexttile;
    plot_grid_with_black_background(koff_heterog(:,:,ads_layer)', c_limits_koff);
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    colorbar;

    % --- Row 3: Plot s_max ---
    smax_slice_heterog = smax_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_smax = [min(smax_slice_heterog(:)), max(smax_slice_heterog(:))];
    if diff(c_limits_smax) < 1e-9; c_limits_smax(2) = c_limits_smax(1) + 1; end

    % Homogeneous s_max
    nexttile;
    plot_grid_with_black_background(smax_homog(:,:,ads_layer)', c_limits_smax);
    xlabel('y-grid index');
    ylabel('s_{max}', 'FontSize', 10);
    
    % Heterogeneous s_max
    nexttile;
    plot_grid_with_black_background(smax_heterog(:,:,ads_layer)', c_limits_smax);
    xlabel('y-grid index');
    set(gca, 'YTickLabel', []);
    colorbar;
    
    % Add a main y-label for the entire layout
    han = gcf();
    han.CurrentAxes = gca();
    ylabel(han.CurrentAxes.Parent, 'x-grid index', 'FontSize',10)
    
    % Save the figure
    disp('Saving parameter fields plot as param_fields.png...');
    print('Adsorption/Testing_kernel/Figures/param_fields.png', '-dpng', '-r300');
end
