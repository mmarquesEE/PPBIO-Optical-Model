function plot_and_save_parameter_fields(kon_homog, koff_homog, smax_homog, ...
                                        kon_heterog, koff_heterog, smax_heterog, ...
                                        ads_x_range, ads_y_range, ads_layer)
    
    % --- Define Publication Style Parameters ---
    target_fig_width_cm = 8.4;
    base_font_size = 8;

    fig = figure('Name', 'Parameter Fields Comparison');
    
    % Use tiledlayout for better control over spacing
    t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Row 1: Plot k_on ---
    kon_slice_heterog = kon_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_kon = [min(kon_slice_heterog(:)), max(kon_slice_heterog(:))];
    if diff(c_limits_kon) < 1e-9; c_limits_kon(2) = c_limits_kon(1) + 1; end
    
    ax1 = nexttile;
    plot_grid_with_black_background(kon_homog(:,:,ads_layer)', c_limits_kon);
    title('Homogeneous', 'FontSize', base_font_size, 'FontWeight', 'normal');
    ylabel('k_{on}', 'FontSize', base_font_size, 'FontWeight', 'bold');
    set(ax1, 'XTickLabel', [], 'FontSize', base_font_size - 1); 
    
    ax2 = nexttile;
    plot_grid_with_black_background(kon_heterog(:,:,ads_layer)', c_limits_kon);
    title('Heterogeneous', 'FontSize', base_font_size, 'FontWeight', 'normal');
    set(ax2, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', base_font_size - 1);
    
    cb1 = colorbar;
    cb1.FontSize = base_font_size - 2;
    cb1.Ticks = linspace(c_limits_kon(1), c_limits_kon(2), 3);
    tick_values_kon = cb1.Ticks;
    scaled_labels_kon = arrayfun(@(x) sprintf('%.1f', x/1000), tick_values_kon, 'UniformOutput', false);
    cb1.TickLabels = scaled_labels_kon;
    cb1.Label.String = '×10^3';
    cb1.Label.FontSize = base_font_size-2;

    % --- Row 2: Plot k_off ---
    koff_slice_heterog = koff_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_koff = [min(koff_slice_heterog(:)), max(koff_slice_heterog(:))];
    if diff(c_limits_koff) < 1e-9; c_limits_koff(2) = c_limits_koff(1) + 1; end
    
    ax3 = nexttile;
    plot_grid_with_black_background(koff_homog(:,:,ads_layer)', c_limits_koff);
    ylabel('k_{off}', 'FontSize', base_font_size, 'FontWeight', 'bold');
    set(ax3, 'XTickLabel', [], 'FontSize', base_font_size - 1); 

    ax4 = nexttile;
    plot_grid_with_black_background(koff_heterog(:,:,ads_layer)', c_limits_koff);
    set(ax4, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', base_font_size - 1);
    
    cb2 = colorbar;
    cb2.FontSize = base_font_size - 2;
    cb2.Ticks = linspace(c_limits_koff(1), c_limits_koff(2), 3);
    tick_values_koff = cb2.Ticks;
    scaled_labels_kon = arrayfun(@(x) sprintf('%.1f', x*1000), tick_values_koff, 'UniformOutput', false);
    cb2.TickLabels = scaled_labels_kon;
    cb2.Label.String = '×10^{-3}';
    cb2.Label.FontSize = base_font_size-2;

    % --- Row 3: Plot s_max ---
    smax_slice_heterog = smax_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    c_limits_smax = [min(smax_slice_heterog(:)), max(smax_slice_heterog(:))];
    if diff(c_limits_smax) < 1e-9; c_limits_smax(2) = c_limits_smax(1) + 1; end
    
    ax5 = nexttile;
    plot_grid_with_black_background(smax_homog(:,:,ads_layer)', c_limits_smax);
    ylabel('s_{max}', 'FontSize', base_font_size, 'FontWeight', 'bold');
    set(ax5, 'FontSize', base_font_size - 1);

    ax6 = nexttile;
    plot_grid_with_black_background(smax_heterog(:,:,ads_layer)', c_limits_smax);
    set(ax6, 'YTickLabel', [], 'FontSize', base_font_size - 1);
    
    cb3 = colorbar;
    cb3.FontSize = base_font_size - 2;
    cb3.Ticks = linspace(c_limits_smax(1), c_limits_smax(2), 4);
    tick_values_smax = cb3.Ticks;
    int_labels = arrayfun(@(x) sprintf('%d', round(x)), tick_values_smax, 'UniformOutput', false);
    cb3.TickLabels = int_labels;

    % --- Add shared labels to the entire layout ---
    xlabel(t, 'y-grid index', 'FontSize', base_font_size);
    ylabel(t, 'x-grid index', 'FontSize', base_font_size);
    
    % --- FIX: Corrected the file path ---
    % Removed leading slash and added missing slash
    save_pub_fig(fig, 'Adsorption/Testing_kernel/Figures/param_fields_pub', target_fig_width_cm);
    close(fig);
end