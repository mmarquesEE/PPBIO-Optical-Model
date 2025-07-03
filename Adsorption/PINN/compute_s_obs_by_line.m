function s_obs_matrix = compute_s_obs_by_line(s_grid, ads_x_range, ads_y_range, ads_layer)
% Computes a matrix of sensorgrams, one for each line in the y-direction.
% Output size: [n_time_points x n_y_lines]
    ads_cells = s_grid(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    % Sum over the x-dimension (dim 2) and the z-dimension (dim 4, which is singleton)
    s_obs_matrix = squeeze(sum(ads_cells, [2, 4]));
end