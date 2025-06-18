function plot_grid_with_black_background(grid_data, c_limits)
    % This function uses transparency to make zero-value areas reveal a black background
    
    % Plot the image and get a handle to it
    h = imagesc(grid_data);
    
    % Set the colormap for the data
    colormap(gca, parula);
    
    % Create a transparency map: 1 for non-zero data, 0 for zero-data
    alpha_map = double(grid_data ~= 0);
    
    % Apply the transparency map
    set(h, 'AlphaData', alpha_map);
    
    % Set the axis background color to black
    set(gca, 'Color', 'k');
    
    % Apply the color limits and tighten the axis
    caxis(c_limits);
    axis tight;
end