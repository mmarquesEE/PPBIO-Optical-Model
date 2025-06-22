function save_pub_fig(fig_handle, file_path_no_ext, target_width_cm)
    % Sets figure properties for publication and saves as PDF and EPS.
    % --- NEW: Automatically creates the output directory if it doesn't exist ---
    
    % Step 1: Extract the folder path from the full file path
    folder_path = fileparts(file_path_no_ext);
    
    % Step 2: If the folder path is not empty and the folder does not exist, create it
    if ~isempty(folder_path) && ~exist(folder_path, 'dir')
        fprintf('Creating output directory: %s\n', folder_path);
        mkdir(folder_path);
    end
    
    % --- The rest of the function is the same as before ---
    set(fig_handle, 'PaperUnits', 'centimeters');
    fig_pos = get(fig_handle, 'Position');
    fig_aspect_ratio = fig_pos(4) / fig_pos(3);
    set(fig_handle, 'Position', [fig_pos(1), fig_pos(2), target_width_cm*40, target_width_cm*40*fig_aspect_ratio]);
    paper_height_cm = target_width_cm * fig_aspect_ratio;
    set(fig_handle, 'PaperSize', [target_width_cm, paper_height_cm]);
    set(fig_handle, 'PaperPosition', [0, 0, target_width_cm, paper_height_cm]);
    
    fprintf('Saving figure to: %s.(pdf/eps)\n', file_path_no_ext);
    print(fig_handle, [file_path_no_ext, '.pdf'], '-dpdf', '-r300');
    print(fig_handle, [file_path_no_ext, '.eps'], '-depsc'); 
end