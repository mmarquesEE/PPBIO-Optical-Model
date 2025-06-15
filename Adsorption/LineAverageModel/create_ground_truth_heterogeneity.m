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