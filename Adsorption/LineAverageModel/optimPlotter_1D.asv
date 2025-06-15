function stop = optimPlotter_1D(log_params, optimValues, state, ...
                               true_params_1D, exp_data, exp_settings, model_config, ads_y_dim)
% optimPlotter_1D visualiza o progresso da otimização com diagnósticos avançados.
    stop = false;
    persistent handles; 
    
    switch state
        case 'init'
            % Na primeira chamada, cria a figura e os eixos
            fig = figure('Name', '1D Optimization Progress & Diagnostics', 'Position', [50, 50, 1800, 900]);
            
            % --- Configuração do Layout 3x4 ---
            % Linha 1: Parâmetros
            handles.ax_kon = subplot(3, 4, 1);
            handles.ax_koff = subplot(3, 4, 2);
            handles.ax_smax = subplot(3, 4, 3);
            
            % Linha 2: Ajuste do Sensorgram
            handles.ax_fit1 = subplot(3, 4, 5);
            handles.ax_fit2 = subplot(3, 4, 6);
            handles.ax_fit3 = subplot(3, 4, 7);
            
            % Linha 3: Resíduos do Sensorgram
            handles.ax_res1 = subplot(3, 4, 9);
            handles.ax_res2 = subplot(3, 4, 10);
            handles.ax_res3 = subplot(3, 4, 11);
            
            % Canto direito: Gráfico Q-Q
            handles.ax_qq = subplot(3, 4, [4, 8, 12]); % Ocupa 3 células na vertical

            % --- Armazena handles e dados ---
            handles.fig = fig;
            handles.true_params_1D = true_params_1D;
            handles.exp_data = exp_data;
            handles.ads_y_dim = ads_y_dim;
            
            % Linka os eixos X do ajuste e do resíduo
            linkaxes([handles.ax_fit1, handles.ax_res1], 'x');
            linkaxes([handles.ax_fit2, handles.ax_res2], 'x');
            linkaxes([handles.ax_fit3, handles.ax_res3], 'x');
            
            updatePlots(log_params, optimValues, handles); % Desenha o estado inicial
        case 'iter'
            if ishandle(handles.fig)
                updatePlots(log_params, optimValues, handles);
            else
                stop = true;
                fprintf('Animation figure closed. Stopping optimization.\n');
            end
        case 'done'
            if ishandle(handles.fig)
                sgtitle(handles.fig, 'Optimization Finished!', 'FontSize', 16, 'FontWeight', 'bold');
            end
    end
end