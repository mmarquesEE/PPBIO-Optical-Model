function LineExperiments()
    clearvars; close all; clc;
    % Shared parameters
    gridN_x = 22; gridN_y = 5; gridN_z = 3;
    ads_layer = 1;
    ads_x_range = [5,15]; ads_y_range = [1,5];

    % --- MODIFIED --- Get number of lines in adsorption region
    ads_y_dim = ads_y_range(2) - ads_y_range(1) + 1;

    % Homogeneous parameters
    homog_params = [9.4e3, 0.0078, 2960]; % kon, koff, smax_total
    
    % Fixed parameters
    D_coeff = 6e-5;
    ru_to_m = 1e-10;
    base_max_velocity = 8.3;
    grid_size_x = 11.0; % mm
    grid_size_z = 0.3; % mm
    
    dx = grid_size_x / gridN_x; % Size of one grid cell in x-direction (mm)
    dz = grid_size_z / gridN_z; % Size of one grid cell in z-direction (mm)
    % Example pulse parameters
    n_exp = 1;
    T1 = 2200; T2 = 2*T1; T3 = 3*T1;
    t_total = 4*T1;
    c_diss = 0;
    c1 = 330e-5;
    c2 = 100e-5;
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
   % =========================================================================
    % --- NOVA SEÇÃO: VISUALIZAÇÃO DAS SUPERFícIES DOS PARÂMETROS ---
    % =========================================================================
    fprintf('\nVisualizando as superfícies de parâmetros geradas...\n');
    
    % Crie uma grade de coordenadas X e Y para o plot 3D
    [X, Y] = meshgrid(1:gridN_y, 1:gridN_x);
    
    % Extraia a camada 2D de cada parâmetro que contém os valores de adsorção
    kon_slice = squeeze(kon_grid_heterog(:, :, ads_layer));
    koff_slice = squeeze(koff_grid_heterog(:, :, ads_layer));
    smax_slice = squeeze(smax_grid_heterog(:, :, ads_layer));
    
    % Crie a figura e os subplots
    figure('Name', 'Visualização dos Parâmetros da Superfície', 'Position', [100, 100, 1600, 500]);
    
    % --- Gráfico de Superfície para k_on ---
    subplot(1, 3, 1);
    surf(X, Y, kon_slice);
    title('Superfície do Parâmetro k_{on}');
    xlabel('Posição da Linha (y)');
    ylabel('Posição ao Longo do Fluxo (x)');
    zlabel('Valor de k_{on}');
    colorbar;
    view(30, 45); % Ajusta o ângulo da câmera para melhor visualização
    
    % --- Gráfico de Superfície para k_off ---
    subplot(1, 3, 2);
    surf(X, Y, koff_slice);
    title('Superfície do Parâmetro k_{off}');
    xlabel('Posição da Linha (y)');
    ylabel('Posição ao Longo do Fluxo (x)');
    zlabel('Valor de k_{off}');
    colorbar;
    view(30, 45);
    
    % --- Gráfico de Superfície para s_max ---
    subplot(1, 3, 3);
    surf(X, Y, smax_slice);
    title('Superfície do Parâmetro s_{max}');
    xlabel('Posição da Linha (y)');
    ylabel('Posição ao Longo do Fluxo (x)');
    zlabel('Valor de s_{max}');
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
    
    % Preallocate matrix to store the calculated resonance angle for each line at each time
    theta_spr_vs_time = zeros(size(s_obs_ru)); formula_response_vs_time = zeros(size(s_obs_ru));
    baseline_offset = calculate_sensorgram_from_formula(n_bulk, n_bulk, n1, d2, wavelength);
    fprintf('Calculated baseline offset of %.4f will be subtracted.\n', baseline_offset);
    tic;
    fprintf('Generating SPR IMAGES FOR VIDEO...\n');
    video_frames_folder = 'spr_video_frames';
    if ~exist(video_frames_folder, 'dir')
        mkdir(video_frames_folder);
    else
        % Opcional: Limpa a pasta de frames antigos antes de começar
        delete(fullfile(video_frames_folder, '*.png'));
    end
    num_time_points = size(s_obs_ru, 1);    % Loop through each time point
    fig_handle = figure('Visible', 'off'); % Figura invisível
    ax = gca;
    ax.Position = [0 0 1 1]; % Eixos preenchem toda a figura
    axis(ax, 'off'); % Remove os eixos e bordas brancas

    % Plota a primeira imagem para inicializar o objeto de imagem e obter um handle
    initial_image_matrix = zeros(ads_y_dim, length(angle_range));
    h_img = imagesc(ax, angle_range, 1:ads_y_dim, initial_image_matrix);
    colormap(ax, 'gray');
    
    % Define os limites de cor uma vez para evitar o auto-ajuste
    % Precisamos calcular o min/max da refletividade em todo o experimento
    % Para simplificar, vamos usar [0, 1], que é o range padrão da refletividade.
    caxis(ax, [0 1]);
    for t_idx = 1:size(s_obs_ru, 1)
        spr_image_matrix = zeros(ads_y_dim, length(angle_range));
        % Loop through each line
        parfor j_idx = 1:size(s_obs_ru, 2)
            % 1. Get the surface concentration in RU for this line at this time
            current_ru = s_obs_ru(t_idx, j_idx);
            
            % 2. Convert RU to the refractive index of the analyte layer (n2)
            n2_analyte = n_bulk + (current_ru*RU_TO_RIU);
            
            % 3. Calculate the full SPR curve and find the resonance angle
            [Rp_curve, resonance_angle] = fresnel_spr_curve(angle_range, n0, n1, n2_analyte, n_bulk, d1, d2, wavelength);
            spr_image_matrix(j_idx, :) = Rp_curve;
            % 4. Store the result
            theta_spr_vs_time(t_idx, j_idx) = resonance_angle;

            % 2. Calculate the sensorgram point using the new function
            delta_neff = calculate_sensorgram_from_formula(n2_analyte, n_bulk, n1, d2, wavelength);
            
            % 3. Store the result
            formula_response_vs_time(t_idx, j_idx) = delta_neff- baseline_offset;
        end
        % Atualize apenas os dados ('CData') do objeto de imagem, não recrie o plot
        set(h_img, 'CData', spr_image_matrix);
        
        % Capture o frame da figura
        frame = getframe(fig_handle);
        
        % Salve o frame usando 'imwrite', que é muito mais rápido que 'saveas'
        filename = fullfile(video_frames_folder, sprintf('frame_%04d.png', t_idx));
        imwrite(frame.cdata, filename);
        % --- Opcional: Exibir um progresso ---
        if mod(t_idx, 100) == 0
            fprintf('Gerado frame %d de %d...\n', t_idx, num_time_points);
        end
    end
    toc;
    close(fig_handle);
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
    fprintf('\n--- Iniciando a criação do vídeo a partir dos frames salvos ---\n');
    tic;

    % --- Configura o objeto VideoWriter ---
    video_filename = 'spri_simulation_final.mp4';
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
end
function s_obs_matrix = compute_s_obs_by_line(s_grid, ads_x_range, ads_y_range, ads_layer)
% Computes a matrix of sensorgrams, one for each line in the y-direction.
% Output size: [n_time_points x n_y_lines]
    ads_cells = s_grid(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    % Sum over the x-dimension (dim 2) and the z-dimension (dim 4, which is singleton)
    s_obs_matrix = squeeze(sum(ads_cells, [2, 4]));
end
function exp_settings = generate_experiments(M, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2)
    % Generates M experiments with orthogonal concentration profiles and flow velocities
    %
    % Inputs:
    %   M - Number of experiments
    %   base_max_velocity - Reference flow velocity (e.g., 8.3)
    %   T1, T2, T3 - Fixed pulse times
    %   t_total - Total experiment duration
    %   c_diss - Dissociation concentration
    %   c1, c2 - Base concentrations
    
    exp_settings = struct(...
        'pulse_times', {}, ...
        'pulse_concs', {}, ...
        't_total', {}, ...
        'max_velocity', {}, ...
        'c_diss', {} ...
    );
    
    % Generate orthogonal concentration pairs using polar coordinates
    angles = linspace(0, pi/2, M);  % Cover quadrant for positive concentrations
    factors = linspace(0.3, 3, M);  % Concentration scaling factors
    
    for i = 1:M
        if M ==1
            M=2;
        end
        % Create orthogonal concentration profiles
        conc_factor1 = factors(ceil(i/2)) * cos(angles(i));
        conc_factor2 = factors(ceil(i/2)) * sin(angles(i));
        
        % Ensure minimum concentration variation
        min_conc = 0.1 * min(c1, c2);
        conc1 = max(c1 * (0.5 + conc_factor1), min_conc);
        conc3 = max(c2 * (0.5 + conc_factor2), min_conc);
        
        % Create velocity profile (logarithmic spacing)
        vel_min = 0.1 * base_max_velocity;
        vel_max = 5.0 * base_max_velocity;
        velocity = exp(log(vel_min) + (i-1)/(M-1) * (log(vel_max) - log(vel_min)));
        
        % Special patterns for every 3rd experiment
        if mod(i,3) == 0
            exp_settings(i).pulse_concs = [conc1, c2, conc3];  % Middle pulse active
        elseif mod(i,4) == 0
            exp_settings(i).pulse_concs = [c1, 0, conc3];      % First pulse fixed
        else
            exp_settings(i).pulse_concs = [conc1, 0, conc3];   % Standard pattern
        end
        
        % Assign common parameters
        exp_settings(i).pulse_times = [T1, T2, T3];
        exp_settings(i).t_total = t_total;
        exp_settings(i).max_velocity = velocity;
        exp_settings(i).c_diss = c_diss;
    end
end
% ================== NEW HELPER FUNCTIONS ==================
% --- Modified compute_residuals_log function ---
function all_residuals = compute_residuals_1D_model(log_params, exp_settings, exp_data, model_config)
% Computes residuals and includes robust error handling for simulation failures.
% This version is optimized to preallocate the results vector for performance.

    p_1D = 10.^log_params;

    % --- START: PREALLOCATION LOGIC ---

    % 1. First, calculate the total number of residual points across all experiments.
    total_num_residuals = 0;
    for exp_idx = 1:length(exp_settings)
        % numel() gets the total number of elements in the data matrix
        total_num_residuals = total_num_residuals + numel(exp_data{exp_idx}.signals);
    end
    
    % 2. Preallocate the full residuals vector with zeros.
    all_residuals = zeros(total_num_residuals, 1);
    
    % 3. Initialize an index to track our position in all_residuals.
    current_idx_start = 1;

    % --- END: PREALLOCATION LOGIC ---


    % --- Part 2: Model-Data Residuals with Error Handling ---
    for exp_idx = 1:length(exp_settings)
        setting = exp_settings(exp_idx);
        s_data_matrix = exp_data{exp_idx}.signals;
        
        try
            [~, s_sim_matrix] = run_single_experiment_1D_model(p_1D, setting, model_config);
            if ~isequal(size(s_sim_matrix), size(s_data_matrix))
                error('Simulation output size mismatch.');
            end
            res_matrix = s_sim_matrix - s_data_matrix;
        catch ME
            fprintf('Warning: Simulation failed. Penalizing parameter set. Error: %s\n', ME.message);
            
            penalty_value = 1e6 * (1 + norm(s_data_matrix(:)));
            res_matrix = penalty_value * ones(size(s_data_matrix));
        end
        
        % --- START: MODIFIED RESULT HANDLING ---
        
        % Calculate the number of elements for this specific experiment
        num_res_in_exp = numel(res_matrix);
        
        % Define the range in the preallocated vector to fill
        idx_range = current_idx_start : (current_idx_start + num_res_in_exp - 1);
        
        % Place the current residuals into the correct slice, ensuring it's a column
        all_residuals(idx_range) = res_matrix(:);
        
        % Update the starting index for the next iteration
        current_idx_start = current_idx_start + num_res_in_exp;
        
        % --- END: MODIFIED RESULT HANDLING ---
    end
end


% ================== EXISTING HELPER FUNCTIONS ==================
function [t, s_obs_by_line] = run_single_experiment_1D_model(p_1D, setting, model_config)
% Runs a simulation using the 1D heterogeneous parameter model.
    ads_y_dim = model_config.ads_y_range(2) - model_config.ads_y_range(1) + 1;
    ads_x_dim = model_config.ads_x_range(2) - model_config.ads_x_range(1) + 1;
    dx = model_config.dx; dz = model_config.dz;
    % Unpack 1D parameter vector
    kon_1D  = p_1D(1:ads_y_dim);
    koff_1D = p_1D(ads_y_dim+1 : 2*ads_y_dim);
    smax_1D = p_1D(2*ads_y_dim+1 : 3*ads_y_dim);

    % Create 2D parameter grids by replicating parameters along the x-axis
    kon_ads  = repmat(kon_1D', ads_x_dim, 1);
    koff_ads = repmat(koff_1D', ads_x_dim, 1);
    smax_ads = repmat(smax_1D', ads_x_dim, 1);
    
    % Create full 3D parameter grids
    [kon_grid, koff_grid, smax_grid] = create_heterogeneous_grids_from_ads(...
        model_config.gridN_x, model_config.gridN_y, model_config.gridN_z, ...
        model_config.ads_x_range, model_config.ads_y_range, model_config.ads_layer, ...
        kon_ads, koff_ads, smax_ads);

    % Create velocity, time, and concentration profiles
    velocity_profile = create_velocity_profile(model_config.gridN_z, setting.max_velocity);
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    s0_grid = zeros(model_config.gridN_x, model_config.gridN_y, model_config.gridN_z);

    % Run simulation
    [t, ~, s] = simulate_3d_flow_model_with_pulses(...
        model_config.gridN_x, model_config.gridN_y, model_config.gridN_z, ...
        kon_grid, koff_grid, smax_grid, velocity_profile, ...
        t_breaks, concentrations, model_config.D_coeff, model_config.ru_to_m, s0_grid,dx,dz);
    
    % Compute the per-line observed signal
    s_obs_by_line = compute_s_obs_by_line(s, model_config.ads_x_range, ...
        model_config.ads_y_range, model_config.ads_layer);
end

function plot_parameter_recovery_1D(p_true, p_opt, p_init, ads_y_dim)
    
    % Extract parameter groups
    true_kon = p_true(1:ads_y_dim);
    true_koff = p_true(ads_y_dim+1:2*ads_y_dim);
    true_smax = p_true(2*ads_y_dim+1:end);
    
    opt_kon = p_opt(1:ads_y_dim);
    opt_koff = p_opt(ads_y_dim+1:2*ads_y_dim);
    opt_smax = p_opt(2*ads_y_dim+1:end);

    init_kon = p_init(1:ads_y_dim);
    init_koff = p_init(ads_y_dim+1:2*ads_y_dim);
    init_smax = p_init(2*ads_y_dim+1:end);
    
    figure('Position', [100, 100, 800, 900]);
    
    % Plot kon recovery
    subplot(3,1,1);
    plot(true_kon, 'ro-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_kon, 'bx--', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_kon, 'g*-', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Line-by-Line Binding Rate (k_{on}) Recovery');
    legend('True (Line Avg)', 'Initial Guess', 'Recovered');
    ylabel('Value'); xlabel('Line Index (j)');
    grid on;

    % Plot koff recovery
    subplot(3,1,2);
    plot(true_koff, 'ro-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_koff, 'bx--', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(opt_koff, 'g*-', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Line-by-Line Unbinding Rate (k_{off}) Recovery');
    ylabel('Value'); xlabel('Line Index (j)');
    grid on;
    
    % Plot smax recovery
    subplot(3,1,3);
    plot(true_smax, 'ro-', 'MarkerSize', 8, 'LineWidth', 2); hold on;
    plot(init_smax, 'bx--', 'MarkerSize', 8, 'LineWidth' ,1.5);
    plot(opt_smax, 'g*-', 'MarkerSize', 8, 'LineWidth', 1.5);
    title('Line-by-Line Maximum Binding (s_{max}) Recovery');
    ylabel('Value'); xlabel('Line Index (j)');
    grid on;
    
    sgtitle('Parameter Recovery Results for 1D Heterogeneous Model');
end

function [kon_grid, koff_grid, smax_grid] = create_heterogeneous_grids_from_ads(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_ads, koff_ads, smax_ads)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_ads;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_ads;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_ads;
end

function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz-1);
    h = nz-1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1,1,nz]);
end

function [kon_grid, koff_grid, smax_grid] = create_ground_truth_heterogeneity(nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    % Parâmetros base
    kon_base = 9.4e3;
    koff_base = 0.0078;
    smax_total = 2960;
    
    % Dimensões da região de adsorção
    num_x_ads = ads_x_range(2) - ads_x_range(1) + 1;
    num_y_ads = ads_y_range(2) - ads_y_range(1) + 1; % Número de linhas

    % Inicializa as matrizes de parâmetros
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    rng(42); % Para resultados aleatórios reprodutíveis

    % --- LÓGICA MODIFICADA PARA HETEROGENEIDADE 1D (POR LINHA) ---

    % 1. Gere um VETOR de valores aleatórios, um para cada LINHA (direção y).
    kon_rand_vec = randn(1, num_y_ads);
    koff_rand_vec = randn(1, num_y_ads);
    smax_rand_vec = randn(1, num_y_ads);

    % 2. Use 'repmat' para replicar o vetor de linha, criando uma matriz
    kon_perturbation = repmat(kon_rand_vec, num_x_ads, 1);
    koff_perturbation = repmat(koff_rand_vec, num_x_ads, 1);
    smax_perturbation = repmat(smax_rand_vec, num_x_ads, 1);

    % 3. Crie as matrizes de parâmetros heterogêneos em 1D
    kon_vals = kon_base * (1 + 0.05 * kon_perturbation);
    koff_vals = koff_base * (1 + 0.05 * koff_perturbation);
    smax_vals = (smax_total / num_y_ads) * (1 + 0.05 * smax_perturbation);
    
    % Atribui os valores à região de adsorção na grade completa
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_vals;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_vals;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_vals;
end

function [t_all, c_s, s] = simulate_3d_flow_model_with_pulses(...
    nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid,dx,dz)
    
    % Initialize state variables
    num_cells = nx * ny * nz;
    c_s0 = zeros(nx, ny, nz);
    c_s0(1, :, :) = concentrations(1); % Initial concentration
    s0 = s0_grid;
    %Q0 = zeros(nx, ny, nz);
    %R0 = zeros(nx, ny, nz);
    y0 = [c_s0(:); s0(:)]; %Q0(:); R0(:)];
    
    % Setup ODE options
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
    
    % --- START: PREALLOCATION LOGIC ---
    
    % 1. Calculate the total number of points for preallocation
    num_segments = length(t_breaks) - 1;
    num_points_per_segment = 1000;
    % Total points = points from segment 1 + points from all other segments (excluding duplicates)
    total_points = num_points_per_segment + (num_segments - 1) * (num_points_per_segment - 1);
    
    % 2. Preallocate results arrays using zeros()
    t_all = zeros(total_points, 1);
    y_all = zeros(total_points, length(y0));
    
    % 3. Initialize an index to keep track of where to insert data
    last_idx = 0;
    
    % --- END: PREALLOCATION LOGIC ---
    
    % Process each time segment
    for seg = 1:num_segments
        t_start = t_breaks(seg);
        t_end = t_breaks(seg+1);
        c0_seg = concentrations(seg);
        
        % Determine time points for segment
        tspan = linspace(t_start, t_end, num_points_per_segment);
        
        % Run simulation for segment
        [t_seg, y_seg] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, ...
            kon_grid, koff_grid, smax_grid, c0_seg, D_coeff, ru_to_m,dx, dz), tspan, y0, options);
        
        % --- START: MODIFIED RESULT HANDLING ---
        
        if seg == 1
            % For the first segment, add all points
            num_to_add = num_points_per_segment;
            current_indices = (last_idx + 1):(last_idx + num_to_add);
            t_all(current_indices) = t_seg;
            y_all(current_indices, :) = y_seg;
            last_idx = last_idx + num_to_add;
        else
            % For subsequent segments, skip the first point to avoid duplicates
            num_to_add = num_points_per_segment - 1;
            current_indices = (last_idx + 1):(last_idx + num_to_add);
            t_all(current_indices) = t_seg(2:end);
            y_all(current_indices, :) = y_seg(2:end, :);
            last_idx = last_idx + num_to_add;
        end

        % --- END: MODIFIED RESULT HANDLING ---
        
        % Update initial condition for next segment
        if seg < num_segments
            y0 = y_seg(end, :)';
            c_s_end = reshape(y0(1:num_cells), [nx, ny, nz]);
            s_end = reshape(y0(num_cells+1:2*num_cells), [nx, ny, nz]);
            c_s_end(1, :, :) = concentrations(seg+1);
            y0 = [c_s_end(:); s_end(:); y0(2*num_cells+1:end)]; 
        end
    end
    
    % Extract variables
    c_s = reshape(y_all(:, 1:num_cells), [length(t_all), nx, ny, nz]);
    s = reshape(y_all(:, num_cells+1:2*num_cells), [length(t_all), nx, ny, nz]);
    %Q = reshape(y_all(:, 2*num_cells+1:3*num_cells), [length(t_all), nx, ny, nz]);
    %R = reshape(y_all(:, 3*num_cells+1:end), [length(t_all), nx, ny, nz]);
    
    % Compute kernel
    %K = exp(-Q) .* R;
end

function dydt = ode_system(~, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0, D_coeff, ru_to_m, dx, dz)
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);
    %Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);
    %R = reshape(y(3*num_cells + 1:4*num_cells), [nx, ny, nz]);
    %dcsdt = zeros(nx, ny, nz);
    %dsdt = zeros(nx, ny, nz);

    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c0;
    dcsdt(1,:,:) = 0;

    % Diffusion terms
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = (c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:))/ (dx^2);
    
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = (c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2))/ (dz^2);
    d2c_dz2(:,:,1) = (c_s(:,:,2) - 2*c_s(:,:,1) + c_s(:,:,1))/ (dz^2);
    d2c_dz2(:,:,end) = (c_s(:,:,end-1) - 2*c_s(:,:,end) + c_s(:,:,end-1))/ (dz^2);
    
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2);
    
    % Advection
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)/ (dx)));
    
    % Adsorption kinetics
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Compute dQ/dt and dR/dt
    %dQdt = kon_grid .* c_s + koff_grid;
    %dRdt = c_s .* exp(Q);

    % Combine all derivatives
    %dydt = [dcsdt(:); dsdt(:); dQdt(:); dRdt(:)];
    dydt = [dcsdt(:); dsdt(:)];
end
% --- NEW FUNCTION ---
function stop = optimPlotter_1D(log_params, optimValues, state, ...
                               true_params_1D, exp_data, exp_settings, model_config, ads_y_dim)
% optimPlotter_1D visualizes the optimization progress for the 1D heterogeneous model.

    stop = false; % This function does not stop the optimization by default
    
    % Use persistent variables to store plot handles
    persistent handles; 
    
    switch state
        case 'init'
            % On the first call, create the figure and axes
            fig = figure('Name', '1D Optimization Progress', 'Position', [100, 100, 1600, 700]);
            
            % --- Setup Axes ---
            ax_kon = subplot(2, 3, 1);
            ax_koff = subplot(2, 3, 2);
            ax_smax = subplot(2, 3, 3);
            
            % We will plot the fit for the first experiment
            ax_sig1 = subplot(2, 3, 4);
            ax_sig2 = subplot(2, 3, 5);
            ax_sig3 = subplot(2, 3, 6);
            
            % --- Store handles and static data ---
            handles.fig = fig;
            handles.ax_kon = ax_kon; handles.ax_koff = ax_koff; handles.ax_smax = ax_smax;
            handles.ax_sigs = [ax_sig1, ax_sig2, ax_sig3];
            
            handles.true_params_1D = true_params_1D;
            handles.exp_data = exp_data;
            handles.exp_settings = exp_settings;
            handles.model_config = model_config;
            handles.ads_y_dim = ads_y_dim;
            
            % Initial plot to set up the view
            updatePlots(log_params, optimValues, handles);

        case 'iter'
            % On each subsequent iteration, update the plots if the figure is still open
            if ishandle(handles.fig)
                updatePlots(log_params, optimValues, handles);
            else
                stop = true;
                fprintf('Animation figure closed. Stopping optimization.\n');
            end
            
        case 'done'
            if ishandle(handles.fig)
                sgtitle(handles.fig, 'Optimization Finished!', 'FontSize', 14, 'FontWeight', 'bold');
            end
    end
end

function updatePlots(current_log_params, optimVals, h)
% Updates plots during optimization, using interpolation to prevent size mismatch errors.

    % --- Update Parameter Plots (this part is unchanged) ---
    current_params_linear = 10.^current_log_params;
    ads_y_dim = h.ads_y_dim;
    % ... (the code for splitting and plotting the parameters kon, koff, smax remains the same) ...
    true_kon = h.true_params_1D(1:ads_y_dim);
    true_koff = h.true_params_1D(ads_y_dim+1 : 2*ads_y_dim);
    true_smax = h.true_params_1D(2*ads_y_dim+1 : end);
    opt_kon = current_params_linear(1:ads_y_dim);
    opt_koff = current_params_linear(ads_y_dim+1 : 2*ads_y_dim);
    opt_smax = current_params_linear(2*ads_y_dim+1 : end);
    cla(h.ax_kon); hold(h.ax_kon, 'on');
    plot(h.ax_kon, true_kon, 'ro-', 'LineWidth', 2, 'DisplayName', 'True');
    plot(h.ax_kon, opt_kon, 'g*-', 'LineWidth', 1.5, 'DisplayName', 'Current');
    title(h.ax_kon, sprintf('k_{on} (Iter: %d)', optimVals.iteration));
    legend(h.ax_kon, 'Location', 'best'); grid(h.ax_kon, 'on'); xlabel(h.ax_kon, 'Line Index');
    cla(h.ax_koff); hold(h.ax_koff, 'on');
    plot(h.ax_koff, true_koff, 'ro-', 'LineWidth', 2);
    plot(h.ax_koff, opt_koff, 'g*-', 'LineWidth', 1.5);
    title(h.ax_koff, sprintf('k_{off} (F-count: %d)', optimVals.funccount));
    grid(h.ax_koff, 'on'); xlabel(h.ax_koff, 'Line Index');
    cla(h.ax_smax); hold(h.ax_smax, 'on');
    plot(h.ax_smax, true_smax, 'ro-', 'LineWidth', 2);
    plot(h.ax_smax, opt_smax, 'g*-', 'LineWidth', 1.5);
    title(h.ax_smax, sprintf('s_{max} (Residual: %.2e)', optimVals.resnorm));
    grid(h.ax_smax, 'on'); xlabel(h.ax_smax, 'Line Index');
    
    % --- MODIFIED: Update Sensorgram Plots with Interpolation ---
    %setting = h.exp_settings(1);
    
    % Unpack the experimental data struct
    data_struct = h.exp_data{1};
    t_exp = data_struct.time;          % The fixed, "master" time vector
    exp_data_matrix = data_struct.signals;

    % Simulate with current parameters to get the new data and its time vector
%     [t_sim, s_sim_matrix] = run_single_experiment_1D_model(current_params_linear, setting, h.model_config);
    % The residual is for ALL experiments, so extract the part for the first one.
    num_points_exp1 = numel(exp_data_matrix);
    residual_exp1 = optimVals.residual(1:num_points_exp1);
    
    % Reconstruct the simulated signal from the residual. It's much faster!
    s_sim_matrix = exp_data_matrix + reshape(residual_exp1, size(exp_data_matrix));

    % Use interpolation to resample the simulated data onto the experimental time grid
    %s_sim_interp = interp1(t_sim, s_sim_matrix, t_exp, 'linear', 'extrap');
    
    % Plot a few sample lines using the common time vector 't_exp'
    lines_to_plot = unique([1, round(ads_y_dim/2), ads_y_dim]);
    for i = 1:length(lines_to_plot)
        line_idx = lines_to_plot(i);
        ax = h.ax_sigs(i);
        
        cla(ax); hold(ax, 'on');
        % Now, both vectors will have the same length (length(t_exp))
        plot(ax, t_exp, exp_data_matrix(:, line_idx), 'b-', 'LineWidth', 2, 'DisplayName', 'Data');
        plot(ax, t_exp, s_sim_matrix(:, line_idx), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Fit');
        title(ax, sprintf('Sensorgram Fit for Line %d', line_idx));
        legend(ax, 'Location', 'best'); xlabel(ax, 'Time (s)'); ylabel(ax, 's_{obs,j}(t)');
        grid(ax, 'on');
    end
    drawnow; % Force the figure window to update
end

% --- NEW HELPER FUNCTION 1: Fresnel Model for 4 Layers ---
function [Rp, resonance_angle] = fresnel_spr_curve(angles_deg, n0, n1, n2, n3, d1, d2, wavelength)
% Calculates the SPR reflectivity curve for p-polarized light in a 4-layer system.
% This is a MATLAB implementation of the provided Python/JAX Fresnel model.

    % Convert input angles from degrees to radians
    th = deg2rad(angles_deg);
    
    % Initialize output array for reflectivity
    Rp = zeros(size(th));
    
    % This loop is the equivalent of jax.vmap, applying the calculation for each angle
    for i = 1:length(th)
        current_th = th(i);
        
        % Check for TIR condition to avoid complex numbers in sqrt where not needed
        sin_th_sq = (n0 * sin(current_th))^2;
        
        % Compute q values for each layer (using complex numbers for generality)
        q0 = sqrt(n0^2 - sin_th_sq + 0i) / n0^2;
        q1 = sqrt(n1^2 - sin_th_sq + 0i) / n1^2;
        q2 = sqrt(n2^2 - sin_th_sq + 0i) / n2^2;
        q3 = sqrt(n3^2 - sin_th_sq + 0i) / n3^2;
        
        % Compute beta values for layers 1 and 2
        beta1 = 2 * pi * d1 * sqrt(n1^2 - sin_th_sq + 0i) / wavelength;
        beta2 = 2 * pi * d2 * sqrt(n2^2 - sin_th_sq + 0i) / wavelength;
        
        % Layer matrices M1 and M2
        M1 = [cos(beta1), -1j * sin(beta1) / q1; 
              -1j * q1 * sin(beta1), cos(beta1)];
          
        M2 = [cos(beta2), -1j * sin(beta2) / q2; 
              -1j * q2 * sin(beta2), cos(beta2)];
        
        % Overall matrix product
        M = M1 * M2;
        
        % Reflection coefficient for p-polarized light
        numerator = (M(1,1) + M(1,2) * q3) * q0 - (M(2,1) + M(2,2) * q3);
        denominator = (M(1,1) + M(1,2) * q3) * q0 + (M(2,1) + M(2,2) * q3);
        rp = numerator / denominator;
        
        % Reflectivity is the squared magnitude of the reflection coefficient
        Rp(i) = abs(rp)^2;
    end
    
    % Find the resonance angle (angle of minimum reflectivity)
    [~, min_idx] = min(Rp);
    resonance_angle = angles_deg(min_idx);
end

% --- NEW HELPER FUNCTION FOR THE ANALYTICAL APPROXIMATION ---
% --- NEW HELPER FUNCTION FOR YOUR SPECIFIED FORMULA ---
function delta_neff = calculate_sensorgram_from_formula(n_analyte, n_bulk, n_metal_complex, d_analyte, lambda)
% Calculates the sensorgram response based on the specific analytical
% formula provided by the user in equations (6) and (7).

    % --- Define terms based on your formula's notation ---
    % Epsilon_2r: Real part of the metal's dielectric constant (layer 2 in a 3-layer model)
    epsilon_2r = real(n_metal_complex^2);
    
    % N3: Refractive index of the analyte layer
    N3_sq = n_analyte^2;
    
    % N4: Refractive index of the bulk/environmental medium
    N4 = n_bulk;
    N4_sq = N4^2;

    % d3: Thickness of the layer causing the change.
    % Based on the physics, we will interpret this as the analyte layer thickness, d2.
    d3 = d_analyte;

    % --- Calculate Equation (6) ---
    term1 = (2 * pi * d3) / lambda;
    
    % Note: (-epsilon_2r * N4^2) will be positive since epsilon_2r for gold is negative.
    term2_numerator = (-epsilon_2r * N4_sq)^(3/2);
    term2_denominator = (epsilon_2r - N4_sq)^2;
    term2 = term2_numerator / term2_denominator;

    term3_numerator = N3_sq - N4_sq;
    term3_denominator = N3_sq; % As per your formula
    term3 = term3_numerator / term3_denominator;

    N_s_eff = term1 * term2 * term3 + N4_sq;
    
    % --- Calculate Equation (7) ---
    delta_neff = N_s_eff - N4;
end