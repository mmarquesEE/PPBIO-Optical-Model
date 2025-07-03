function save_subplot_as_eps(ax_handle, file_path)
    % Create a new invisible figure
    temp_fig = figure('Visible', 'off', 'PaperPositionMode', 'auto');
    % Copy the target subplot to the new figure
    new_ax = copyobj(ax_handle, temp_fig);
    % Adjust the subplot's position to fill the new figure
    set(new_ax, 'Position', get(groot, 'defaultAxesPosition'));
    % Save the new figure containing only the subplot
    print(temp_fig, file_path, '-depsc');
    % Close the temporary figure
    close(temp_fig);
end