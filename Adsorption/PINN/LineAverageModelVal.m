function LineAverageModelVal()
    clearvars; close all; clc;
    %% Shared parameters
    generate_video_frames = true; % Mude para 'false' para pular a criação do vídeo e acelerar o script
    gridN_x = 22; gridN_y = 5; gridN_z = 3;ads_layer = 1;
    ads_x_range = [5,15]; ads_y_range = [1,5];
    % Get number of lines in adsorption region
    ads_y_dim = ads_y_range(2) - ads_y_range(1) + 1;
    % Homogeneous parameters
    homog_params = [9.4e3, 0.0078, 2960]; % kon, koff, smax_total
    % Fixed parameters
    D_coeff = 6e-5;ru_to_m = 1e-10;base_max_velocity = 8.3;
    grid_size_x = 11.0; % mm
    grid_size_z = 0.3; % mm
    dx = grid_size_x / gridN_x; % Size of one grid cell in x-direction (mm)
    dz = grid_size_z / gridN_z; % Size of one grid cell in z-direction (mm)
    % Example pulse parameters
    n_exp = 1;T1 = 2200; T2 = 2*T1; T3 = 3*T1;t_total = 4*T1;
    c_diss = 0;c1 = 330e-4;c2 = 100e-4;
    exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);
    model_config.gridN_x = gridN_x;
    model_config.gridN_y = gridN_y;
    model_config.gridN_z = gridN_z;
    model_config.dx = dx; model_config.dz = dz;
    model_config.ads_x_range = ads_x_range;
    model_config.ads_y_range = ads_y_range;
    model_config.ads_layer = ads_layer;
    model_config.D_coeff = D_coeff;
    model_config.ru_to_m = ru_to_m;
    % Create ground truth 2D heterogeneous parameters
    [kon_grid_heterog, koff_grid_heterog, smax_grid_heterog] = ...
        create_ground_truth_heterogeneity(gridN_x, gridN_y, gridN_z, ads_x_range, ads_y_range, ads_layer);
    %=========================================================================
    %% --- New Section: Visualizing the surfaces ---
    % =========================================================================
    fprintf('\nVisualizing the surfaces...\n');
    
    % Crie uma grade de coordenadas X e Y para o plot 3D
    [X, Y] = meshgrid(1:gridN_y, 1:gridN_x);
    
    % Extraia a camada 2D de cada parâmetro que contém os valores de adsorção
    kon_slice = squeeze(kon_grid_heterog(:, :, ads_layer));
    koff_slice = squeeze(koff_grid_heterog(:, :, ads_layer));
    smax_slice = squeeze(smax_grid_heterog(:, :, ads_layer));
    
    % Crie a figura e os subplots
    figure('Name', 'Surface Parameters', 'Position', [100, 100, 1600, 500]);
    
    % --- Gráfico de Superfície para k_on ---
    subplot(1, 3, 1);
    surf(X, Y, kon_slice);
    title('Surface for k_{on}');
    xlabel('Line Position (y)');
    ylabel('Line Position (x)');
    zlabel('k_{on}');
    colorbar;
    view(30, 45); % Ajusta o ângulo da câmera para melhor visualização
    
    % --- Gráfico de Superfície para k_off ---
    subplot(1, 3, 2);
    surf(X, Y, koff_slice);
    title('Surface for k_{off}');
    xlabel('Line Position (y)');
    ylabel('Line Position (x)');
    zlabel('k_{off}');
    colorbar;
    view(30, 45);
    
    % --- Gráfico de Superfície para s_max ---
    subplot(1, 3, 3);
    surf(X, Y, smax_slice);
    title('Surface for s_{max}');
    xlabel('Line Position (y)');
    ylabel('Line Position (x)');
    zlabel('s_{max}');
    colorbar;
    view(30, 45);
    
    sgtitle('Visualização 3D da Heterogeneidade dos Parâmetros na Superfície', 'FontSize', 16, 'FontWeight', 'bold');
    
    % =========================================================================
    % --- STEP 1: GENERATE "EXPERIMENTAL" DATA FIRST ---
    % By creating this data upfront, we know the exact size of all outputs
    % and can reuse the clean data later.
    % =========================================================================
    fprintf('\nGenerating ground-truth data for all experiments...\n');
    tic;
    exp_data = cell(n_exp, 1);
    total_rows = 0;
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        velocity_profile = create_velocity_profile(gridN_z, setting.max_velocity);
        t_breaks = [0, setting.pulse_times, setting.t_total];
        concentrations = [setting.pulse_concs, setting.c_diss];
        s0_grid = zeros(gridN_x, gridN_y, gridN_z);
        
        [t_exp, ~, s_heterog_exp] = simulate_3d_flow_model_with_pulses(...
            gridN_x, gridN_y, gridN_z, kon_grid_heterog, koff_grid_heterog, smax_grid_heterog, ...
            velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid,dx,dz);
        
        s_obs_by_line_clean = compute_s_obs_by_line(s_heterog_exp, ads_x_range, ads_y_range, ads_layer);
        
        noise_level = 0.01;
        noise_matrix = 1 + noise_level * randn(size(s_obs_by_line_clean));
        
        data_struct.time = t_exp;
        data_struct.signals_clean = s_obs_by_line_clean; % Store the clean version
        data_struct.signals = s_obs_by_line_clean .* noise_matrix; % And the noisy version
        exp_data{exp_idx} = data_struct;
        %Plotting
        t_plot = exp_data{exp_idx}.time;
        s_plot_by_line = exp_data{exp_idx}.signals_clean; % Plot the clean data
        s_plot_global = sum(s_plot_by_line, 2);
        
        figure('Position', [100, 100, 1200, 600]);
        subplot(1,2,1);
        plot(t_plot, s_plot_global, 'r-', 'LineWidth', 2);
        title('Global (Summed) Sensorgram'); xlabel('Time (s)'); ylabel('Total s_{obs}(t)'); grid on;
        
        subplot(1,2,2);
        plot(t_plot, s_plot_by_line, 'LineWidth', 1.5);
        title('Line-by-Line Sensorgrams (1D Heterogeneity)'); xlabel('Time (s)'); ylabel('s_{obs, j}(t)');
        legend(arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), 'Location', 'best'); grid on;
        sgtitle(sprintf('Experiment %d: Ground Truth Data', exp_idx));
        % 1. Calculate total rows by summing up data points from all experiments
        total_rows = total_rows + numel(exp_data{exp_idx}.signals_clean);
    end
    fprintf('Data generation complete.\n');    
    toc;
    % ============== IDENTIFIABILITY ANALYSIS (1D HETEROGENEITY) ================
    fprintf('\nStarting 1D Identifiability Analysis...\n');
    
    % --- MODIFIED --- Create the "true" 1D parameter vector
    % We assume the true line parameter is the average over the flow direction
    kon_ads_2D = kon_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    koff_ads_2D = koff_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    smax_ads_2D = smax_grid_heterog(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);

    % Average parameters along the x-axis (flow direction) to get line parameters
    p_true_kon_1D = mean(kon_ads_2D, 1)';   % [ads_y_dim x 1]
    p_true_koff_1D = mean(koff_ads_2D, 1)'; % [ads_y_dim x 1]
    p_true_smax_1D = mean(smax_ads_2D, 1)'; % [ads_y_dim x 1]
    
    p_true_1D = [p_true_kon_1D; p_true_koff_1D; p_true_smax_1D];
    N_params_1D = length(p_true_1D);
    
    % --- START: PREALLOCATION FOR JACOBIAN ---
    tic;
    % 2. Preallocate the combined Jacobian matrix
    J_combined = zeros(total_rows, N_params_1D);
    
    % 3. Initialize a row indexer
    current_row_start = 1;
    
    % --- END: PREALLOCATION FOR JACOBIAN ---
    
    % Loop to compute and fill the Jacobian
    for exp_idx = 1:n_exp
        setting = exp_settings(exp_idx);
        fprintf('Computing Jacobian for Experiment %d (1D model)...\n', exp_idx);
        
        % REUSE baseline output from the data we already generated. No need to re-run.
        s_obs_matrix = exp_data{exp_idx}.signals_clean;
        s_obs_vector = s_obs_matrix(:); % Vectorize for Jacobian
        
        n_rows_exp = length(s_obs_vector);
        J_exp = zeros(n_rows_exp, N_params_1D);
        h = 1e-5;
        
        parfor k = 1:N_params_1D
            p_pert = p_true_1D;
            p_pert(k) = p_pert(k) * (1 + h);
            
            % We only need to run the perturbed simulation here
            [~, s_pert_matrix] = run_single_experiment_1D_model(p_pert, setting, model_config);
            s_pert_vector = s_pert_matrix(:);
            
            J_exp(:, k) = (s_pert_vector - s_obs_vector) / (p_true_1D(k) * h);
        end
        
        % --- Fill the preallocated matrix ---
        row_range = current_row_start : (current_row_start + n_rows_exp - 1);
        J_combined(row_range, :) = J_exp;
        current_row_start = current_row_start + n_rows_exp;
    end
    toc;
    % ================= START OF SVD ANALYSIS =================

    rankJ = rank(J_combined);
    fprintf('\n1D Identifiability Analysis Results:\n');
    fprintf('Total Parameters (3 * N_y): %d\n', N_params_1D);
    fprintf('Rank of Combined Jacobian: %d\n', rankJ);
    
    if rankJ < N_params_1D
        fprintf('WARNING: The model is structurally unidentifiable. Rank < Number of Parameters.\n');
    else
        fprintf('SUCCESS: The model appears to be structurally identifiable (Jacobian has full rank).\n');
    end
    fprintf('Now performing SVD analysis to investigate practical identifiability...\n');
    
    % --- Step 1: Perform Singular Value Decomposition ---
    % Use the 'econ' flag for efficiency, as we only need the first N_params_1D vectors
    [~, S, V] = svd(J_combined, 'econ');
    
    % Extract the diagonal singular values
    singular_values = diag(S);
    
    % --- Step 2: Analyze and Plot Singular Values ---
    figure('Name', 'SVD Analysis of Jacobian', 'Position', [100, 100, 1400, 600]);
    
    subplot(1, 2, 1);
    semilogy(singular_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    title('Singular Values of the Jacobian');
    xlabel('Singular Value Index');
    ylabel('Magnitude (log scale)');
    xlim([0, N_params_1D + 1]);
    % Add text for the condition number
    cond_number = singular_values(1) / singular_values(end);
    legend(sprintf('Condition Number: %.2e', cond_number));
    
    % --- Step 3: Analyze and Plot Parameter Combinations (Right Singular Vectors) ---
    subplot(1, 2, 2);
    imagesc(abs(V)); % Use absolute value for clarity of magnitude
    colorbar;
    title('Parameter Contributions to Singular Vectors (V)');
    xlabel('Singular Vector Index (1=Most Identifiable -> N=Least Identifiable)');
    ylabel('Parameter Index');
    
    % Create meaningful labels for the y-axis
    param_labels = [arrayfun(@(i) sprintf('kon_{%d}', i), 1:ads_y_dim, 'UniformOutput', false), ...
                    arrayfun(@(i) sprintf('koff_{%d}', i), 1:ads_y_dim, 'UniformOutput', false), ...
                    arrayfun(@(i) sprintf('smax_{%d}', i), 1:ads_y_dim, 'UniformOutput', false)];
    yticks(1:N_params_1D);
    yticklabels(param_labels);
    
    sgtitle('SVD-based Identifiability Analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    % ================= END OF SVD ANALYSIS =================

   % ============== PARAMETER IDENTIFICATION (1D HETEROGENEITY) ================    
    fprintf('\nStarting 1D Parameter Identification...\n');
    
    % --- MODIFIED --- Setup initial guess and bounds for 1D model
    homog_kon = homog_params(1);
    homog_koff = homog_params(2);
    homog_smax_per_line = (homog_params(3) / ads_y_dim); % Evenly distributed total smax
    
    p0_base_kon = log10(homog_kon) * ones(ads_y_dim, 1);
    p0_base_koff = log10(homog_koff) * ones(ads_y_dim, 1);
    p0_base_smax = log10(homog_smax_per_line) * ones(ads_y_dim, 1);
    
    p0_base = [p0_base_kon; p0_base_koff; p0_base_smax];
    
        % Define TIGHT bounds in log10 space to keep the optimizer in a stable region.
    % These correspond to kon=[1e2, 1e7], koff=[1e-5, 1], smax_per_line=[1e-4, 5]
    lb_kon = log10(1e3);  ub_kon = log10(1e7);
    lb_koff = log10(1e-6); ub_koff = log10(1e-1); 
    lb_smax = log10(1e-4); ub_smax = log10(3000);   
    
    lb = [repmat(lb_kon, ads_y_dim, 1); ...
          repmat(lb_koff, ads_y_dim, 1); ...
          repmat(lb_smax, ads_y_dim, 1)];
    ub = [repmat(ub_kon, ads_y_dim, 1); ...
          repmat(ub_koff, ads_y_dim, 1); ...
          repmat(ub_smax, ads_y_dim, 1)];
    % 3. Create the final p0 by adding a small random perturbation AND clamping to the bounds
    rng('default'); % For reproducible randomness
    noise_level = 0.1; % Small perturbation (e.g., 0.1 standard deviations in log space)
    p0_rand = p0_base + noise_level * randn(size(p0_base));
    
    % Clamp the randomized p0 to be within the bounds
    p0 = max(lb, p0_rand); % Enforce lower bound
    p0 = min(ub, p0);     % Enforce upper bound
    % --- NEW --- Create a handle to the plotter function with all necessary data
    optim_plot_fun = @(log_p, optim_v, state) optimPlotter_1D(...
        log_p, optim_v, state, ...
        p_true_1D, exp_data, exp_settings, model_config, ads_y_dim);

    % --- MODIFIED --- Add the 'OutputFcn' to your optimization options
    optim_opts = optimoptions('lsqnonlin', ...
        'Algorithm', 'trust-region-reflective', ...
        'Display', 'iter', ...
        'MaxIterations', 50, ...
        'UseParallel', true, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10, ...
        'OutputFcn', optim_plot_fun); % This tells lsqnonlin to call our plotter
        
    % --- MODIFIED --- The residual function remains the same
    residual_fun = @(log_params) compute_residuals_1D_model(log_params, exp_settings, exp_data, model_config);
    tic;
    % Run optimization
    [opt_log_params, ~] = lsqnonlin(residual_fun, p0, lb, ub, optim_opts);
    toc;
    % --- MODIFIED --- Analyze and plot results for 1D model
    opt_params_1D = 10.^opt_log_params;
    
    % Plot recovery of the 1D parameters
    plot_parameter_recovery_1D(p_true_1D, opt_params_1D, 10.^p0, ads_y_dim);
    % =========================================================================
    % --- FINAL, PHYSICALLY-ACCURATE VALIDATION WORKFLOW ---
    % =========================================================================
    fprintf('\n--- Starting Full Physical Model Validation ---\n');

    % --- Step 1: Define Optical and Physical Constants ---
    fprintf('Defining optical parameters...\n');
    wavelength = 670; % nm
    d1 = 50;          % Gold film thickness (nm)
    % The analyte layer (n2) thickness is effectively infinite for the evanescent wave
    d2 = 1000;        % Effectively infinite analyte layer (nm)
    
    n0 = sqrt(2.3104);         % Optical substrate (Prism)
    n1 = sqrt(-14.379 + 1.0084j); % Gold film (complex RI)
    n_bulk = sqrt(1.7876);         % Flow cell solution (baseline buffer)
    
    % Define the angular range for SPR curve calculation
    angle_range = linspace(65, 80, 500); % [start_angle, end_angle, num_points]
    
    % Define the conversion factor from Response Units (RU) to Refractive Index Units (RIU)
    % 1000 RU = 0.001 RIU change
    RU_TO_RIU = 0.001 / 1000;
    
    % --- Step 2: Get the Ground-Truth Sensorgram Data (in RU) ---
    % We use the clean, line-averaged data from our initial simulation
    t_exp = exp_data{1}.time;
    s_obs_ru = exp_data{1}.signals_clean; % Sensorgrams in RU

    % --- Step 3: Convert Sensorgrams to Resonance Angles via Fresnel Model ---
    fprintf('Processing %d time points for %d lines...\n', size(s_obs_ru, 1), size(s_obs_ru, 2));
    
    % Preallocate matrices to store the calculated results
    theta_spr_vs_time = zeros(size(s_obs_ru));
    formula_response_vs_time = zeros(size(s_obs_ru));
    
    baseline_offset = calculate_sensorgram_from_formula(n_bulk, n_bulk, n1, d2, wavelength);
    fprintf('Calculated baseline offset of %.4f will be subtracted.\n', baseline_offset);
    video_frames_folder = 'spr_video_frames';
    % --- Setup for Video Frame Generation (only if flag is true) ---
    if generate_video_frames
        fprintf('Generating SPR IMAGES FOR VIDEO...\n');
        if ~exist(video_frames_folder, 'dir')
            mkdir(video_frames_folder);
        else
            delete(fullfile(video_frames_folder, '*.png'));
        end
        
        % Create the hidden figure and axes ONCE before the loop
        fig_handle = figure('Visible', 'off');
        ax = gca;
        ax.Position = [0 0 1 1];
        axis(ax, 'off');
        initial_image_matrix = zeros(ads_y_dim, length(angle_range));
        h_img = imagesc(ax, angle_range, 1:ads_y_dim, initial_image_matrix);
        colormap(ax, 'gray');
        caxis(ax, [0 1]);
    end
    
    % --- Main Calculation and Optional Frame Generation Loop ---
    tic;
    num_time_points = size(s_obs_ru, 1);
    fprintf('Processing %d time points...\n', num_time_points);
    n2_vs_time = zeros(size(s_obs_ru)); 

    for t_idx = 1:size(s_obs_ru, 1)
        spr_image_matrix = zeros(ads_y_dim, length(angle_range));
        % Loop through each line
       parfor j_idx = 1:size(s_obs_ru, 2)
            current_ru = s_obs_ru(t_idx, j_idx);
            n2_analyte = n_bulk + (current_ru * RU_TO_RIU);
            n2_vs_time(t_idx, j_idx) = n2_analyte; % Salva o RI calculado
            % Calculate full SPR curve, resonance angle, and formula response
            [Rp_curve, resonance_angle] = fresnel_spr_curve(angle_range, n0, n1, n2_analyte, n_bulk, d1, d2, wavelength);
            spr_image_matrix(j_idx, :) = Rp_curve;
            theta_spr_vs_time(t_idx, j_idx) = resonance_angle;
            delta_neff = calculate_sensorgram_from_formula(n2_analyte, n_bulk, n1, d2, wavelength);
            formula_response_vs_time(t_idx, j_idx) = delta_neff - baseline_offset;
        end
        % --- This block for creating and saving images ONLY runs if the flag is true ---
        if generate_video_frames
            % Update the data of the existing image object
            set(h_img, 'CData', spr_image_matrix);
            
            % Capture the frame
            frame = getframe(fig_handle);
            
            % Save the frame efficiently
            filename = fullfile(video_frames_folder, sprintf('frame_%04d.png', t_idx));
            imwrite(frame.cdata, filename);
        end
        
        % Display progress
        if mod(t_idx, 100) == 0
            fprintf('Processed frame %d de %d...\n', t_idx, num_time_points);
        end
    end
    % --- Clean up the figure and print completion message (only if flag was true) ---
    if generate_video_frames
        close(fig_handle);
        fprintf('All %d frames were generated and saved.\n', num_time_points);
    end
    toc;
    figure('Name', 'SPR Curves for All Lines at Max Response');
    hold on;
    for j_idx = 1:ads_y_dim
        plot(angle_range, spr_image_matrix(j_idx, :), 'LineWidth', 2, 'DisplayName', sprintf('Linha %d', j_idx));
    end
    hold off;
    grid on;
    xlabel('Ângulo de Incidência (graus)');
    ylabel('Refletividade');
    legend('Location', 'best');
    ylim([0, 1]);

     % =========================================================================
    % --- ETAPA 2: CRIAR O VÍDEO A PARTIR DOS FRAMES SALVOS ---
    % =========================================================================
    if generate_video_frames
        fprintf('\n--- Iniciando a criação do vídeo a partir dos frames salvos ---\n');
        tic;
    
        % --- Configura o objeto VideoWriter ---
        video_filename = 'Adsorption/LineAverageModel/spri_simulation_final.mp4';
        outputVideo = VideoWriter(video_filename, 'MPEG-4');
        outputVideo.FrameRate = 30;
        open(outputVideo);
    
        % --- Pega a lista de todos os arquivos de imagem ---
        image_files_struct = dir(fullfile(video_frames_folder, '*.png'));
        image_files_cell = {image_files_struct.name};
        
        % --- Ordena os nomes dos arquivos numericamente (robusto) ---
        % Extrai os números dos nomes dos arquivos
        str_nums = regexp(image_files_cell, '\d+', 'match', 'once');
        num_vals = str2double(str_nums);
        
        % Ordena os números e pega os índices da ordenação
        [~, sorted_indices] = sort(num_vals);
        
        % Usa os índices para ordenar a lista de nomes de arquivos
        sorted_image_files = image_files_cell(sorted_indices);
        
        % --- Loop através dos arquivos ordenados para escrever o vídeo ---
        fprintf('Lendo %d frames para criar o vídeo...\n', length(sorted_image_files));
        for i = 1:length(sorted_image_files)
            % Monta o caminho completo para o arquivo de imagem
            img_path = fullfile(video_frames_folder, sorted_image_files{i});
            
            % Lê a imagem
            img = imread(img_path);
            
            % Escreve o frame no vídeo
            writeVideo(outputVideo, img);
        end
    
        % --- Finaliza e fecha o arquivo de vídeo ---
        close(outputVideo);
        toc;
        fprintf('\nVídeo salvo com sucesso como "%s".\n', video_filename);
    end

    % --- Step 4: Plot the "Proper" Sensorgram and Validate ---
    % We plot the Resonance Angle directly. To validate, we overlay
    % the original RU data on a second y-axis to show the shapes match.
    
    figure('Name', 'Proper Sensorgram: Resonance Angle vs. Time', 'Position', [300, 300, 1400, 700]);
    
    lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);
    
    for i = 1:length(lines_to_plot)
        subplot(1, length(lines_to_plot), i);
        line_idx = lines_to_plot(i);
        
        % --- This is the proper, physically-correct sensorgram ---
        plot(t_exp, theta_spr_vs_time(:, line_idx), 'r-', 'LineWidth', 2, 'DisplayName', 'Sensorgram (Resonance Angle)');
        
        grid on;
        xlabel('Time (s)');
        ylabel('Resonance Angle (degrees)');
        title(sprintf('Proper Sensorgram for Line %d', line_idx));
        
        % --- For validation, plot the original RU data on a separate y-axis ---
        yyaxis right % Activate the right y-axis
        plot(t_exp, s_obs_ru(:, line_idx), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Original Simulation (RU)');
        ylabel('Response Units (RU)');
        
        legend('Location', 'best');
        xlim([0, t_exp(end)]);
    end
    sgtitle('Final Validation: The Physically Correct Sensorgram (Angle vs. Time)', 'FontSize', 16);

    % --- Step 4: Plot the Sensorgram from Your Formula ---
    % =========================================================================
    % --- FINAL PLOTTING OF RESULTS ---
    % =========================================================================
    
    % --- Step 1: Calculate the Global Sensorgram ---
    % We sum the responses from all lines at each time point by summing along
    % the second dimension (the columns) of the results matrix.
    global_sensorgram = sum(formula_response_vs_time, 2);

    % --- Step 2: Create the Figure and Subplots ---
    figure('Name', 'Final Sensorgram Results from Formula', 'Position', [200, 200, 1400, 600]);

    % --- Plot 1: All Individual Line Sensorgrams ---
    subplot(1, 2, 1);
    plot(t_exp, formula_response_vs_time, 'LineWidth', 1.5);
    grid on;
    title('Individual Line Sensorgrams');
    xlabel('Time (s)');
    ylabel('Change from Baseline (\Delta{N}_s^{eff})');
    xlim([0, t_exp(end)]);
    % Optional: Add a legend if you have a small number of lines
    if ads_y_dim <= 10
        legend(arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), 'Location', 'best');
    end

    % --- Plot 2: Global (Summed) Sensorgram ---
    subplot(1, 2, 2);
    plot(t_exp, global_sensorgram, 'r-', 'LineWidth', 2);
    grid on;
    title('Global (Summed) Sensorgram');
    xlabel('Time (s)');
    ylabel('Total Change (\Sigma \Delta{N}_s^{eff})');
    xlim([0, t_exp(end)]);
    
    sgtitle('Final Sensorgrams Calculated from Analytical Formula', 'FontSize', 16, 'FontWeight', 'bold');
    
    absolute_neff_vs_time = formula_response_vs_time + n_bulk;

    % --- Step 2: Calculate the Global (Summed) Absolute N_s_eff ---
    global_absolute_neff = n_bulk+global_sensorgram;

    % --- Step 3: Create the Figure and Subplots ---
    figure('Name', 'Absolute Effective RI Sensorgrams', 'Position', [200, 200, 1400, 600]);

    % --- Plot 1: All Individual Line Sensorgrams ---
    subplot(1, 2, 1);
    plot(t_exp, absolute_neff_vs_time, 'LineWidth', 1.5);
    grid on;
    title('Individual Line Sensorgrams (Absolute N_s^{eff})');
    xlabel('Time (s)');
    ylabel('Absolute Effective RI ({N}_s^{eff})');
    xlim([0, t_exp(end)]);
    % Optional: Add a legend if you have a small number of lines
    if ads_y_dim <= 10
        legend(arrayfun(@(j) sprintf('Line %d', j), 1:ads_y_dim, 'UniformOutput', false), 'Location', 'best');
    end

    % --- Plot 2: Global (Summed) Sensorgram ---
    subplot(1, 2, 2);
    plot(t_exp, global_absolute_neff, 'r-', 'LineWidth', 2);
    grid on;
    title('Global (Summed) Sensorgram');
    xlabel('Time (s)');
    ylabel('Absolute Effective RI (\Sigma {N}_s^{eff})');
    xlim([0, t_exp(end)]);
    
    sgtitle('Final Sensorgrams Plotted as Absolute N_s^{eff}', 'FontSize', 16, 'FontWeight', 'bold');
    % =========================================================================
    % --- VALIDAÇÃO FINAL: O PROCESSO INVERSO ---
    % =========================================================================
    
    % --- Step 1: Analisa os frames salvos para reconstruir o sensorgram ---
    % Esta função executa o processo inverso que você descreveu
    theta_extracted_from_frames = analyze_spr_frames_to_get_sensorgram(video_frames_folder, angle_range);
    % --- Step 2: Converte os Ângulos Extraídos de volta para Índice de Refração ---
    fprintf('Convertendo ângulos extraídos de volta para Índice de Refração via interpolação...\n');
    
    n2_reconstructed_vs_time = zeros(size(theta_extracted_from_frames));
    
    % Para cada linha, use a relação (n2 vs theta) que calculamos como uma tabela de consulta
    parfor j_idx = 1:ads_y_dim
        % Os dados conhecidos da nossa simulação "forward"
        known_n2s = n2_vs_time(:, j_idx);
        known_thetas = theta_spr_vs_time(:, j_idx);
        
        % Os ângulos que queremos converter de volta para n2
        angles_to_convert = theta_extracted_from_frames(:, j_idx);

        % A relação theta vs n2 é monotônica, mas para robustez, removemos pontos duplicados
        [unique_thetas, unique_indices] = unique(known_thetas);
        unique_n2s = known_n2s(unique_indices);

        % Use interp1 para fazer a "inversão" da função
        n2_reconstructed_vs_time(:, j_idx) = interp1(unique_thetas, unique_n2s, angles_to_convert, 'linear', 'extrap');
    end
    % --- Step 3: Crie o gráfico de comparação final ---
    % Este é o teste definitivo: os dados originais vs. os dados re-extraídos.
    % Ambos estão em unidades de "graus", então a comparação é direta.
    
    figure('Name', 'Validação Final: Sensorgram Original vs. Extraído dos Frames', 'Position', [300, 300, 1400, 700]);
    
    lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);

    for i = 1:length(lines_to_plot)
        subplot(1, length(lines_to_plot), i);
        line_idx = lines_to_plot(i);
        
        % Plota os ângulos de ressonância calculados DIRETAMENTE da física
        plot(t_exp, theta_spr_vs_time(:, line_idx), 'b-', 'LineWidth', 4, 'DisplayName', 'Original (Cálculo Direto)');
        hold on;
        
        % Plota os ângulos de ressonância EXTRAÍDOS da análise das imagens salvas
        plot(t_exp, theta_extracted_from_frames(:, line_idx), 'r--', 'LineWidth', 2, 'DisplayName', 'Extraído das Imagens');
        
        grid on;
        xlabel('Time (s)');
        ylabel('Resonance Angle (degrees)');
        title(sprintf('Validação Final para a Linha %d', line_idx));
        legend('Location', 'best');
        xlim([0, t_exp(end)]);
    end
    figure('Name', 'Validação Final: RI Original vs. RI Reconstruído', 'Position', [300, 300, 1400, 700]);

    for i = 1:length(lines_to_plot)
            subplot(1, length(lines_to_plot), i);
            line_idx = lines_to_plot(i);
            
            % Plota o perfil de RI original (calculado a partir de s_obs_ru)
            plot(t_exp, n2_vs_time(:, line_idx), 'b-', 'LineWidth', 4, 'DisplayName', 'RI Original (da Simulação)');
            hold on;
            
            % Plota o perfil de RI reconstruído a partir da análise das imagens
            plot(t_exp, n2_reconstructed_vs_time(:, line_idx), 'r--', 'LineWidth', 2, 'DisplayName', 'RI Reconstruído (das Imagens)');
            
            grid on;
            xlabel('Time (s)');
            ylabel('Refractive Index (RIU)');
            title(sprintf('Validação de RI para a Linha %d', line_idx));
            legend('Location', 'best');
            xlim([0, t_exp(end)]);
   end
    sgtitle('Validação Completa do Modelo e do Processo de Análise', 'FontSize', 16);

end 