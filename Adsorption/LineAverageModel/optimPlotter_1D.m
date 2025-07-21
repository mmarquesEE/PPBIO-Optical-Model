function stop = optimPlotter_1D(log_params, optimValues, state, true_params_1D, exp_data, ads_y_dim)
    % This plotter function visualizes the optimization progress and saves it as an MP4 video.
    
    % Use persistent variables to hold handles and data across iterations
    persistent fig_handle plot_data video_writer; 
    
    switch state
        case 'init'
            % --- Initialization Step (runs once at the beginning) ---
            
            % 1. Store static data that won't change during the optimization
            plot_data.true_params_1D = true_params_1D;
            plot_data.exp_data = exp_data;
            plot_data.ads_y_dim = ads_y_dim;
            
            % 2. Call updatePlots to create and draw the initial figure
            fig_handle = updatePlots(log_params, optimValues, plot_data);
            
            % 3. Set up the VideoWriter object to save the animation
            video_filename = 'Adsorption/LineAverageModel/Videos/optimization_progress.mp4';
            if ~exist(fileparts(video_filename), 'dir')
               mkdir(fileparts(video_filename));
            end
            video_writer = VideoWriter(video_filename, 'MPEG-4');
            video_writer.FrameRate = 2; % Adjust for desired speed (e.g., 2-5 fps)
            video_writer.Quality = 100; % Set quality (0-100, 100 is best)
            
            % 4. Open the video file for writing
            open(video_writer);
            
            % 5. Capture and write the very first frame (initial guess)
            frame = getframe(fig_handle);
            writeVideo(video_writer, frame);
            
        case 'iter'
            % --- Iteration Step (runs after each algorithm step) ---
            
            % 1. Check if the figure window was closed by the user
            if ~ishandle(fig_handle)
                stop = true; % Stop the optimizer
                close(video_writer); % Ensure the video file is closed
                return;
            end
            
            % 2. Refresh the plots with the current parameter values
            updatePlots(log_params, optimValues, plot_data);
            
            % 3. Capture the updated figure as a new frame
            frame = getframe(fig_handle);
            
            % 4. Write the new frame to the video file
            writeVideo(video_writer, frame);
            
        case 'done'
            % --- Finalization Step (runs once at the end) ---
            
            % 1. Check if the figure and video objects still exist
            if ishandle(fig_handle) && ~isempty(video_writer)
                % Save the final figure as a static image for reference
                % save_pub_fig(fig_handle, 'Adsorption/LineAverageModel/Figures/figure_fit_diagnostics_FINAL', 18);
                
                % 2. Close the video file to finalize and save it
                close(video_writer);
                fprintf('\nOptimization video successfully saved to: %s\n', video_writer.Filename);
            end
    end
    
    stop = false; % By default, tell the optimizer to continue
end