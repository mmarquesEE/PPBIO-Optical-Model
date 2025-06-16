function theta_from_frames = analyze_spr_frames_to_get_sensorgram(frames_folder, angle_range)
% Reads a folder of SPR image frames, finds the resonance angle for each
% line in each frame, and returns the resulting sensorgram matrix.

    fprintf('\n--- Iniciando o processo inverso: Analisando frames para extrair o sensorgram ---\n');
    
    % --- Pega a lista de todos os arquivos de imagem ---
    image_files_struct = dir(fullfile(frames_folder, '*.png'));
    if isempty(image_files_struct)
        error('A pasta de frames "%s" está vazia ou não foi encontrada. Por favor, execute o script com "generate_video_frames = true" primeiro.', frames_folder);
    end
    image_files_cell = {image_files_struct.name};
    
    % --- Ordena os nomes dos arquivos numericamente para garantir a ordem correta ---
    str_nums = regexp(image_files_cell, '\d+', 'match', 'once');
    num_vals = str2double(str_nums);
    [~, sorted_indices] = sort(num_vals);
    sorted_image_files = image_files_cell(sorted_indices);
    
    % --- Lê a primeira imagem para determinar as dimensões ---
    first_img_path = fullfile(frames_folder, sorted_image_files{1});
    first_img = imread(first_img_path);
%     [num_lines, ~] = size(rgb2gray(first_img)); % Get size from grayscale version
    [num_lines, ~] = size(im2gray(first_img)); 
    num_frames = length(sorted_image_files);
    
    % --- Pré-aloca a matriz de resultados para o sensorgram ---
    theta_from_frames = zeros(num_frames, num_lines);
    
    fprintf('Analisando %d frames para %d linhas...\n', num_frames, num_lines);
    tic;
    
    % --- Loop através de cada arquivo de imagem ordenado ---
    for i = 1:num_frames
        img_path = fullfile(frames_folder, sorted_image_files{i});
        img_uint8 = imread(img_path); % imread irá carregar a imagem como uint8
        
        img_gray = im2gray(img_uint8);
        
        % O resto da função continua como antes
        reflectivity_matrix = double(img_gray) / 255.0;
        % --- MUDANÇA AQUI: Loop através de cada linha para aplicar a interpolação ---
%         for j = 1:num_lines
%             rp_curve = reflectivity_matrix(j, :);
%             theta_from_frames(i, j) = find_subpixel_minimum(rp_curve, angle_range);
%         end
        [~, min_indices] = min(reflectivity_matrix, [], 2);
        resonance_angles_for_frame = angle_range(min_indices);
        theta_from_frames(i, :) = resonance_angles_for_frame;
    end
    toc;
    fprintf('Extração do sensorgram a partir das imagens concluída.\n');
end