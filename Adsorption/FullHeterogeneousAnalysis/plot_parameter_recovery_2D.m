function plot_parameter_recovery_2D(p_true, p_opt, p_init, model_config)
    % This function visualizes the 2D parameter recovery using heatmaps.
    
    % --- De-vectorize all parameter sets ---
    num_sites = (model_config.ads_x_range(2) - model_config.ads_x_range(1) + 1) * ...
                (model_config.ads_y_range(2) - model_config.ads_y_range(1) + 1);
    
    true_kon = reshape(p_true(1:num_sites), [model_config.ads_x_range(2)-model_config.ads_x_range(1)+1, model_config.ads_y_range(2)-model_config.ads_y_range(1)+1]);
    opt_kon  = reshape(p_opt(1:num_sites), size(true_kon));
    init_kon = reshape(p_init(1:num_sites), size(true_kon));

    true_koff = reshape(p_true(num_sites+1:2*num_sites), size(true_kon));
    opt_koff  = reshape(p_opt(num_sites+1:2*num_sites), size(true_kon));
    
    true_smax = reshape(p_true(2*num_sites+1:end), size(true_kon));
    opt_smax  = reshape(p_opt(2*num_sites+1:end), size(true_kon));

    % --- Create Figure ---
    fig = figure('Name', '2D Parameter Recovery Results');
    t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- Plot k_on recovery ---
    c_lims_kon = [min(true_kon(:)), max(true_kon(:))];
    nexttile; imagesc(true_kon'); caxis(c_lims_kon); axis xy; title('True k_{on}'); ylabel('k_{on}');
    nexttile; imagesc(opt_kon');  caxis(c_lims_kon); axis xy; title('Estimated k_{on}');
    nexttile; imagesc(abs(true_kon - opt_kon)'); axis xy; title('Absolute Error'); colorbar;
    
    % --- Plot k_off recovery ---
    c_lims_koff = [min(true_koff(:)), max(true_koff(:))];
    nexttile; imagesc(true_koff'); caxis(c_lims_koff); axis xy; title('True k_{off}'); ylabel('k_{off}');
    nexttile; imagesc(opt_koff');  caxis(c_lims_koff); axis xy; title('Estimated k_{off}');
    nexttile; imagesc(abs(true_koff - opt_koff)'); axis xy; title('Absolute Error'); colorbar;
    
    % --- Plot s_max recovery ---
    c_lims_smax = [min(true_smax(:)), max(true_smax(:))];
    nexttile; imagesc(true_smax'); caxis(c_lims_smax); axis xy; title('True s_{max}'); ylabel('s_{max}');
    nexttile; imagesc(opt_smax');  caxis(c_lims_smax); axis xy; title('Estimated s_{max}');
    nexttile; imagesc(abs(true_smax - opt_smax)'); axis xy; title('Absolute Error'); colorbar;
    
    save_pub_fig(fig, 'Adsorption/FullHeterogeneousAnalysis/Figures/figure_param_recovery_2D', 12);
    close(fig);
end