function save_pub_fig(fig_handle, file_path_no_ext, target_width_cm)
    % Sets figure properties for publication and saves as PDF and EPS.
    set(fig_handle, 'PaperUnits', 'centimeters');
    pos = get(fig_handle, 'Position');
    aspect_ratio = pos(4) / pos(3);
    target_height_cm = target_width_cm * aspect_ratio;
    set(fig_handle, 'PaperSize', [target_width_cm, target_height_cm]);
    set(fig_handle, 'PaperPosition', [0, 0, target_width_cm, target_height_cm]);
    print(fig_handle, [file_path_no_ext, '.pdf'], '-dpdf', '-r300');
    print(fig_handle, [file_path_no_ext, '.eps'], '-depsc');
    set(fig_handle, 'Color', 'none'); 
    set(fig_handle, 'InvertHardcopy', 'off');
    print(fig_handle, [file_path_no_ext, '.svg'], '-dsvg');
    fprintf('Saved figure to %s.pdf and %s.eps\n', file_path_no_ext, file_path_no_ext);
end