function transport_limited_adsorption
    clear all;close all;clc;
    % Parameters 
    N = 2;                          % Number of binding sites
    koff = [1e-3, 0.015];           % Dissociation rates (s^-1)
    KD = [1e-9, 30e-9];             % Equilibrium constants (M)
    kon = koff ./ KD;               % Calculate association rate constants
    smax = [100, 100];              % Maximum capacities (RU)
    k_tr = 1e8;                     % Transport rate (RU/M/s)
    c0_values = [1, 3, 10, 30, 100, 200] * 1e-9; % nM -> M
    
    % Time span 
    t_assoc = 1000;                 % Association phase duration (s)
    t_diss = 1000;                  % Dissociation phase duration (s)
    
    % Initialize figure handles
    figure; ax1 = gca; hold on; title('Site 1 Binding Traces');
    figure; ax2 = gca; hold on; title('Site 2 Binding Traces');
    figure; ax3 = gca; hold on; title('Total Binding Traces');
    figure; ax4 = gca; hold on; title('Surface Compartment Concentration');

    % Colors for different concentrations
    colors = lines(length(c0_values));
    
    % Simulate for each concentration
    for i = 1:length(c0_values)
        % ===== Association Phase =====
        y0_assoc = [0; zeros(N, 1)]; % Initial conditions
        [t1, y1] = ode15s(@(t, y) ode_system(t, y, kon, koff, smax, k_tr, c0_values(i)), ...
                          [0 t_assoc], y0_assoc);
        
        % ===== Dissociation Phase =====
        y0_diss = y1(end, :)';       % Final state of association
        [t2, y2] = ode15s(@(t, y) ode_system(t, y, kon, koff, smax, k_tr, 0), ...
                          [0 t_diss], y0_diss);
        
        % Combine results and add noise
        c_s = [y1(:, 1); y2(2:end, 1)]; % Surface compartment concentration
        t_full = [t1; t1(end) + t2(2:end)];  % Exclude t=0 from dissociation
        s1 = [y1(:, 2); y2(2:end, 2)] + randn(length(t_full), 1)*1; 
        s2 = [y1(:, 3); y2(2:end, 3)] + randn(length(t_full), 1)*1;
        total_s = s1 + s2;
        
        % Plot results
        % Plotting
        plot(ax1, t_full, s1, 'Color', colors(i, :), 'LineWidth', 1.5);
        plot(ax2, t_full, s2, 'Color', colors(i, :), 'LineWidth', 1.5);
        plot(ax3, t_full, total_s, 'Color', colors(i, :), 'LineWidth', 1.5);
        plot(ax4, t_full, c_s, 'LineWidth', 1.5);
    end
    
    % Format plots
    format_plot(ax1, 'Site 1 Coverage (RU)');
    format_plot(ax2, 'Site 2 Coverage (RU)');
    format_plot(ax3, 'Total Coverage (RU)');
    format_plot(ax4, 'c_s (M)');
    legend(ax3, arrayfun(@(c) sprintf('%.0f nM', c*1e9), c0_values, 'UniformOutput', false));
    legend(ax4, arrayfun(@(c) sprintf('%.0f nM', c*1e9), c0_values, 'UniformOutput', false));
end
function dydt = ode_system(~, y, kon, koff, smax, k_tr, c0)
    c_s = y(1);
    s = y(2:end);
    
    % Calculate ds/dt for each site (units: RU/s)
    dsdt = kon' .* c_s .* (smax' - s) - koff' .* s;
    
    % Calculate dc_s/dt (units: M/s using conversion factor 1e-6)
    sum_dsdt = sum(dsdt);
    dcsdt = (k_tr * (c0 - c_s) - sum_dsdt) * 1e-6; % Scale for unit consistency
    
    % Combine derivatives
    dydt = [dcsdt; dsdt];
end

function format_plot(ax, ylabel_text)
    grid(ax, 'on');
    xlabel(ax, 'Time (s)');
    ylabel(ax, ylabel_text);
    set(ax, 'FontSize', 12);
end