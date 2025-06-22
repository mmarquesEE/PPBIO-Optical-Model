function save_pub_fig(fig_handle, file_path_no_ext, target_width_cm)
    % Sets figure properties for publication and saves as PDF and EPS.
    set(fig_handle, 'PaperUnits', 'centimeters');
    fig_pos = get(fig_handle, 'Position');
    fig_aspect_ratio = fig_pos(4) / fig_pos(3);
    set(fig_handle, 'Position', [fig_pos(1), fig_pos(2), target_width_cm*40, target_width_cm*40*fig_aspect_ratio]);
    paper_height_cm = target_width_cm * fig_aspect_ratio;
    set(fig_handle, 'PaperSize', [target_width_cm, paper_height_cm]);
    set(fig_handle, 'PaperPosition', [0, 0, target_width_cm, paper_height_cm]);
    print(fig_handle, [file_path_no_ext, '.pdf'], '-dpdf', '-r300');
    print(fig_handle, [file_path_no_ext, '.eps'], '-depsc'); 
    fprintf('Saved publication figure to %s.pdf and .eps\n', file_path_no_ext);
end