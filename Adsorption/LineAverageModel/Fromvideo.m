%% SCRIPT FOR INVERSE PROCESSING OF EXPERIMENTAL SPRI VIDEO
% MODIFIED TO INCLUDE TEMPORAL SMOOTHING TO REDUCE GLITCHES

clearvars; close all; clc;

% =========================================================================
%% --- STEP 1: USER CONFIGURATION & PHYSICAL CONSTANTS ---
% =========================================================================
% (This section remains unchanged)
video_filename = 'Adsorption/LineAverageModel/Videos/spri_simulation_final.avi';
wavelength = 670; d1 = 50; d2 = 1000;
n0 = sqrt(2.3104);
n1 = sqrt(-14.379 + 1.0084j);
n_bulk = sqrt(1.7876);
RU_TO_RIU = 0.001 / 1000;
angle_range = linspace(65, 80, 1000);
fprintf('Configuration loaded for video: %s\n', video_filename);

% =========================================================================
%% --- STEP 2: READ VIDEO AND EXTRACT FRAMES ---
% =========================================================================
% (This section remains unchanged)
fprintf('\n--- Step 2: Reading video file and extracting frames... ---\n');
tic;
try
    videoObj = VideoReader(video_filename);
    fprintf('✅ VideoReader successfully created the object.\n');
catch ME
    fprintf(2, '\n--- DETAILED MATLAB ERROR INFORMATION ---\n');
    fprintf(2, 'Error Message:     %s\n', ME.message);
    fprintf(2, 'Error Identifier:  %s\n', ME.identifier);
    if ~isempty(ME.stack), fprintf(2, 'Error occurred in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line); end
    fprintf(2, '--- END OF ERROR --- \n\n');
    error('Execution stopped because the video file could not be read.');
end
numFrames = videoObj.NumFrames;
vidHeight = videoObj.Height;
vidWidth = videoObj.Width;
frameRate = videoObj.FrameRate;
if vidWidth ~= length(angle_range), error('Video width (%d) does not match angle_range length (%d).', vidWidth, length(angle_range)); end
fprintf('Video details: %d frames, %d lines (height), %d angle points (width).\n', numFrames, vidHeight, vidWidth);
all_frames = zeros(vidHeight, vidWidth, numFrames, 'uint8');
frame_count = 0;
while hasFrame(videoObj)
    frame_count = frame_count + 1;
    frame_rgb = readFrame(videoObj);
    if size(frame_rgb, 3) == 3, all_frames(:, :, frame_count) = rgb2gray(frame_rgb); else, all_frames(:, :, frame_count) = frame_rgb; end
    if mod(frame_count, 100) == 0, fprintf('Processed frame %d of %d...\n', frame_count, numFrames); end
end
fprintf('Successfully read %d frames.\n', frame_count);
toc;

% =========================================================================
%% --- STEP 3: EXTRACT RESONANCE ANGLE FROM EACH FRAME ---
% =========================================================================
% (This section remains unchanged - it generates the noisy angle data)
fprintf('\n--- Step 3: Extracting resonance angles from frames... ---\n');
tic;
theta_vs_time = zeros(numFrames, vidHeight);
for t_idx = 1:numFrames
    current_frame = all_frames(:, :, t_idx);
    [~, min_indices] = min(current_frame, [], 2);
    theta_vs_time(t_idx, :) = angle_range(min_indices)';
end
fprintf('Resonance angle extraction complete.\n');
toc;

%% <<< NEW SECTION TO SMOOTH THE DATA >>>
% =========================================================================
%% --- STEP 3.5: APPLY TEMPORAL SMOOTHING TO REDUCE JUMPS ---
% =========================================================================
fprintf('\n--- Step 3.5: Applying temporal smoothing filter... ---\n');
tic;

% Define the window size for the moving average filter.
% A larger window = more smoothing. Start with 3 or 5.
smoothing_window_size = 3; 

% Use the smoothdata function to apply the filter to the angle data.
% This will average each point with its neighbors, reducing single-frame glitches.
theta_smoothed = smoothdata(theta_vs_time, 'movmean', smoothing_window_size);

fprintf('Smoothing complete using a %d-frame window.\n', smoothing_window_size);
toc;
%% <<< END OF NEW SECTION >>>

% =========================================================================
%% --- STEP 4: GENERATE CALIBRATION CURVE (ANGLE vs. REFRACTIVE INDEX) ---
% =========================================================================
% (This section remains unchanged)
fprintf('\n--- Step 4: Generating Angle-to-RI calibration curve... ---\n');
tic;
n2_calibration_range = linspace(n_bulk, n_bulk + 0.001, 1000);
theta_calibration = zeros(size(n2_calibration_range));
% Assuming 'fresnel_spr_curve.m' exists on your MATLAB path
parfor i = 1:length(n2_calibration_range)
    n2_current = n2_calibration_range(i);
    [~, resonance_angle] = fresnel_spr_curve(angle_range, n0, n1, n2_current, n_bulk, d1, d2, wavelength);
    theta_calibration(i) = resonance_angle;
end
[theta_calibration_unique, unique_idx] = unique(theta_calibration);
n2_calibration_unique = n2_calibration_range(unique_idx);
fprintf('Calibration curve generation complete.\n');
toc;

% =========================================================================
%% --- STEP 5: INVERSE CONVERSION (ANGLE -> RI -> RU) ---
% =========================================================================
fprintf('\n--- Step 5: Converting angles to RU via interpolation... ---\n');
tic;

%% <<< MODIFIED LINE: Use the smoothed angle data for conversion >>>
% We now use 'theta_smoothed' instead of the original 'theta_vs_time'
n2_reconstructed_vs_time = interp1(theta_calibration_unique, n2_calibration_unique, theta_smoothed, 'linear', 'extrap');

ru_reconstructed = (n2_reconstructed_vs_time - n_bulk) / RU_TO_RIU;
fprintf('Inverse conversion complete.\n');
toc;

% =========================================================================
%% --- STEP 6: PLOT RECONSTRUCTED SENSORGRAMS ---
% =========================================================================
fprintf('\n--- Step 6: Plotting results... ---\n');
time_axis = (0:numFrames-1) / frameRate;
time_axis = linspace(0, 8800, numFrames); % Creates a time axis from 0 to 8800 seconds

figure('Name', 'Smoothed Reconstructed Sensorgrams', 'Position', [100, 100, 900, 600]);
plot(time_axis, ru_reconstructed, 'LineWidth', 1.5);
grid on; box on;

%% <<< MODIFIED LINE: Updated title to reflect smoothing >>>
title('Smoothed Reconstructed Sensorgrams from Experimental Video');
xlabel('Time (s)');
ylabel('Response (RU)');
xlim([0, time_axis(end)]);
legend_labels = arrayfun(@(j) sprintf('Line %d', j), 1:vidHeight, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northwest');
fprintf('Plotting of smoothed data complete.\n');

%% <<< NEW SECTION: OPTIONAL COMPARISON PLOT >>>
% =========================================================================
%% --- STEP 7: COMPARISON OF ORIGINAL VS. SMOOTHED DATA ---
% =========================================================================
% This plot shows the effect of the smoothing filter on a single line
figure('Name', 'Smoothing Filter Comparison', 'Position', [200, 200, 900, 600]);
line_to_plot = 3; % Choose which line to display (e.g., Line 3)

plot(time_axis, theta_vs_time(:, line_to_plot), 'Color', [0.5 0.5 1.0], 'DisplayName', 'Original (Noisy) Angle');
hold on;
plot(time_axis, theta_smoothed(:, line_to_plot), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Angle');
hold off;

grid on; box on;
title(sprintf('Effect of Smoothing Filter on Line %d', line_to_plot));
xlabel('Time (s)');
ylabel('Resonance Angle (degrees)');
legend('Location', 'best');