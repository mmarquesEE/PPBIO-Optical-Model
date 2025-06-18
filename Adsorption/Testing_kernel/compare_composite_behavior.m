function compare_composite_behavior(s_homog, s_heterog, Q_homog, Q_heterog, K_homog, K_heterog,...
        t_homog, t_heterog, ads_x_range, ads_y_range, ads_layer,...
        kon_homog, smax_per_cell, kon_heterog_ads, smax_heterog_ads)
    
    % Extract adsorption region for homogeneous
    ads_cells_s_homog = s_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_Q_homog = Q_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_K_homog = K_homog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    % Extract adsorption region for heterogeneous
    ads_cells_s_heterog = s_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_Q_heterog = Q_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    ads_cells_K_heterog = K_heterog(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    
    s0_homog_ads = squeeze(ads_cells_s_homog(1, :, :));
    s0_heterog_ads = squeeze(ads_cells_s_heterog(1, :, :));
    
    alpha_homog_ads = kon_homog * smax_per_cell * ones(size(kon_heterog_ads));
    alpha_heterog_ads = kon_heterog_ads .* smax_heterog_ads;
    
    alpha_homog_ads = reshape(alpha_homog_ads, [1, size(alpha_homog_ads)]);
    alpha_heterog_ads = reshape(alpha_heterog_ads, [1, size(alpha_heterog_ads)]);
    s0_homog_ads = reshape(s0_homog_ads, [1, size(s0_homog_ads)]);
    s0_heterog_ads = reshape(s0_heterog_ads, [1, size(s0_heterog_ads)]);

    decay_term_homog = s0_homog_ads .* exp(-ads_cells_Q_homog);
    decay_term_heterog = s0_heterog_ads .* exp(-ads_cells_Q_heterog);
    
    kernel_term_homog = alpha_homog_ads .* ads_cells_K_homog;
    kernel_term_heterog = alpha_heterog_ads .* ads_cells_K_heterog;
    
    decay_sum_homog = squeeze(sum(decay_term_homog, [2,3]));
    kernel_sum_homog = squeeze(sum(kernel_term_homog, [2,3]));
    s_obs_homog = decay_sum_homog + kernel_sum_homog;
    
    decay_sum_heterog = squeeze(sum(decay_term_heterog, [2,3]));
    kernel_sum_heterog = squeeze(sum(kernel_term_heterog, [2,3]));
    s_obs_heterog = decay_sum_heterog + kernel_sum_heterog;
    
    discrepancy = trapz(t_heterog, (s_obs_homog - s_obs_heterog).^2);
    fprintf('Composite behavior discrepancy: %.2e\n', discrepancy);
    
    % Plot decomposition for both cases
    fig1 = figure;
    subplot(2,1,1);
    plot(t_homog, s_obs_homog, 'k-', 'LineWidth', 2); hold on;
    plot(t_homog, decay_sum_homog, 'b--', 'LineWidth', 1.5);
    plot(t_homog, kernel_sum_homog, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('s_{obs}');
    legend('Total', 'Decay Term', 'Kernel Term', 'Location', 'best');
    %title('(A) Homogeneous Case: Signal Decomposition');
    title('(A)')
    grid on;
    
    subplot(2,1,2);
    plot(t_heterog, s_obs_heterog, 'k-', 'LineWidth', 2); hold on;
    plot(t_heterog, decay_sum_heterog, 'b--', 'LineWidth', 1.5);
    plot(t_heterog, kernel_sum_heterog, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('s_{obs}');
    legend('Total', 'Decay Term', 'Kernel Term', 'Location', 'best');
    %title('(B) Heterogeneous Case: Signal Decomposition');
    title('(B)')
    grid on;
    
    disp('Saving signal decomposition plot as signal_decomposition.png...');
    print(fig1, 'Adsorption/Testing_kernel/Figures/signal_decomposition.png', '-dpng', '-r300');
    
    % Plot total signal comparison
    fig2 = figure;
    plot(t_heterog, s_obs_homog, 'b-', 'LineWidth', 2); hold on;
    plot(t_heterog, s_obs_heterog, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('s_{obs}(t)');
    legend('Homogeneous', 'Heterogeneous (5% var)', 'Location', 'best');
    %title(sprintf('Composite Behavior Comparison\nDiscrepancy: %.2e', discrepancy));
    grid on;
    
    disp('Saving sensorgram comparison plot as sensorgram_comparison.png...');
    print(fig2, 'Adsorption/Testing_kernel/Figures/sensorgram_comparison.png', '-dpng', '-r300');
end