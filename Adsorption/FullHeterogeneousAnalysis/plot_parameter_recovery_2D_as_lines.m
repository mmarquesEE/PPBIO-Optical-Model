function plot_parameter_recovery_2D_as_lines(p_true, p_opt, p_init, num_sites)
    site_indices = 1:num_sites;
    true_kon = p_true(1:num_sites); opt_kon = p_opt(1:num_sites);
    true_koff = p_true(num_sites+1:2*num_sites); opt_koff = p_opt(num_sites+1:2*num_sites);
    true_smax = p_true(2*num_sites+1:end); opt_smax = p_opt(2*num_sites+1:end);
    
    fig = figure('Name', '2D Parameter Recovery (Line Plots)');
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    nexttile; plot(site_indices, true_kon, 'ro-'); hold on; plot(site_indices, opt_kon, 'k.-'); hold off; grid on; title('k_{on} Recovery'); ylabel('Value'); legend({'True', 'Estimated'});
    nexttile; plot(site_indices, true_koff, 'ro-'); hold on; plot(site_indices, opt_koff, 'k.-'); hold off; grid on; title('k_{off} Recovery'); ylabel('Value');
    nexttile; plot(site_indices, true_smax, 'ro-'); hold on; plot(site_indices, opt_smax, 'k.-'); hold off; grid on; title('s_{max} Recovery'); ylabel('Value'); xlabel('Site Index');
    
    save_pub_fig(fig, 'Adsorption/FullHeterogeneousAnalysis/Figures/figure_param_recovery_2D_lines', 8.4);
    close(fig);
end