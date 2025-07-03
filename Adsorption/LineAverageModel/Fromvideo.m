%% SCRIPT FOR INVERSE PROCESSING OF EXPERIMENTAL SPRI VIDEO
%
% This script reads an experimental SPRI video, extracts the resonance angle
% from each frame for each line, and performs an inverse calculation to
% reconstruct the sensorgrams in Response Units (RU).
%
% Created based on the inverse process outlined in the user's 7-step workflow.

clearvars; close all; clc;

% =========================================================================
%% --- STEP 1: USER CONFIGURATION & PHYSICAL CONSTANTS ---
% =========================================================================
% ⚠️ ACTION REQUIRED: Update the parameters in this section to match your
% experimental setup.

% --- Input Video File ---
% Replace this with the path to your video file.
video_filename = 'Adsorption/LineAverageModel/Videos/received_data_original.avi';

% --- Optical & Physical Parameters ---
% These MUST match the system used to generate the video data.
% These values are taken from Step 4 of your original script.
wavelength = 670; % Wavelength in nanometers
d1 = 50;          % Gold film thickness in nanometers
n0 = sqrt(2.3104);         % Prism refractive index
n1 = sqrt(-14.379 + 1.0084j); % Gold film complex RI
n_bulk = sqrt(1.7876);     % Flow cell solution (baseline buffer) RI
RU_TO_RIU = 0.001 / 1000;  % Conversion factor (1000 RU = 0.001 RIU change)

% --- Angle Range (CRITICAL) ---
% This MUST match the angle range represented by the pixel columns of your video.
% The number of points (last argument) must equal the video's width in pixels.
% Example: If your video is 1000 pixels wide, covering 65 to 80 degrees:
angle_range = linspace(65, 80, 1000); % [start_angle, end_angle, num_pixels]

fprintf('Configuration loaded for video: %s\n', video_filename);

% =========================================================================
%% --- STEP 2: READ VIDEO AND EXTRACT FRAMES ---
% =========================================================================
fprintf('\n--- Step 2: Reading video file and extracting frames... ---\n');
tic;

% Create a VideoReader object
try
    videoObj = VideoReader(video_filename);
catch ME
    error('Failed to open video file: %s. Please check the path and file format.', video_filename);
end

% Get video properties
numFrames = videoObj.NumFrames;
vidHeight = videoObj.Height; % This is the number of sensor lines
vidWidth = videoObj.Width;
frameRate = videoObj.FrameRate;

% --- Sanity Check ---
if vidWidth ~= length(angle_range)
    error('Video width (%d pixels) does not match the length of angle_range (%d points). Please correct the angle_range definition in Step 1.', vidWidth, length(angle_range));
end

fprintf('Video details: %d frames, %d lines (height), %d angle points (width).\n', numFrames, vidHeight, vidWidth);

% Pre-allocate a 3D matrix to store all frames for performance
all_frames = zeros(vidHeight, vidWidth, numFrames, 'uint8');
frame_count = 0;

while hasFrame(videoObj)
    frame_count = frame_count + 1;
    frame_rgb = readFrame(videoObj);
    
    % Convert frame to grayscale if it is in color
    if size(frame_rgb, 3) == 3
        all_frames(:, :, frame_count) = rgb2gray(frame_rgb);
    else
        all_frames(:, :, frame_count) = frame_rgb;
    end
    
    if mod(frame_count, 100) == 0
        fprintf('Processed frame %d of %d...\n', frame_count, numFrames);
    end
end
fprintf('Successfully read %d frames.\n', frame_count);
toc;

% =========================================================================
%% --- STEP 3: EXTRACT RESONANCE ANGLE FROM EACH FRAME ---
% =========================================================================
% This step finds the position of the SPR dip for each line over time.
fprintf('\n--- Step 3: Extracting resonance angles from frames... ---\n');
tic;

% Pre-allocate matrix to store the resonance angle for each frame and line
% Rows = time (frames), Columns = lines
theta_vs_time = zeros(numFrames, vidHeight);

for t_idx = 1:numFrames
    % Get the current 2D frame (lines x angles)
    current_frame = all_frames(:, :, t_idx);
    
    % Find the index of the minimum intensity for all lines at once.
    % The '2' dimension indicates to find the min across each row (columns).
    [~, min_indices] = min(current_frame, [], 2);
    
    % Map these pixel indices to the corresponding angles from angle_range
    % and store them in the correct time-slice (row) of the output matrix.
    % Note the transpose (') to make it a row vector.
    theta_vs_time(t_idx, :) = angle_range(min_indices)';
end

fprintf('Resonance angle extraction complete.\n');
toc;

% =========================================================================
%% --- STEP 4: GENERATE CALIBRATION CURVE (ANGLE vs. REFRACTIVE INDEX) ---
% =========================================================================
% To convert the measured angles back to RI, we need a model of the relationship.
% We generate this by simulating the SPR curve for a range of known RIs.
fprintf('\n--- Step 4: Generating Angle-to-RI calibration curve... ---\n');
tic;

% Define a high-resolution range of analyte RIs to simulate
% This range should comfortably bracket all expected values.
n2_calibration_range = linspace(n_bulk, n_bulk + 0.01, 500);
theta_calibration = zeros(size(n2_calibration_range));

% Use parfor for potential speed-up if you have Parallel Computing Toolbox
parfor i = 1:length(n2_calibration_range)
    n2_current = n2_calibration_range(i);
    
    % For each RI, calculate the full SPR curve and find its resonance angle
    [~, resonance_angle] = calculate_fresnel_spr(angle_range, n0, n1, n2_current, d1, wavelength);
    theta_calibration(i) = resonance_angle;
end

% --- Remove any non-unique points to ensure a valid interpolation table ---
[theta_calibration_unique, unique_idx] = unique(theta_calibration);
n2_calibration_unique = n2_calibration_range(unique_idx);

fprintf('Calibration curve generation complete.\n');
toc;

% =========================================================================
%% --- STEP 5: INVERSE CONVERSION (ANGLE -> RI -> RU) ---
% =========================================================================
fprintf('\n--- Step 5: Converting angles to RU via interpolation... ---\n');
tic;

% Use the calibration curve to convert the measured angles into RIU.
% interp1 looks up each angle in `theta_vs_time` within the calibration
% table and returns the corresponding refractive index.
n2_reconstructed_vs_time = interp1(theta_calibration_unique, n2_calibration_unique, theta_vs_time, 'linear', 'extrap');

% Convert the reconstructed refractive index back to Response Units (RU)
ru_reconstructed = (n2_reconstructed_vs_time - n_bulk) / RU_TO_RIU;

fprintf('Inverse conversion complete.\n');
toc;


% =========================================================================
%% --- STEP 6: PLOT RECONSTRUCTED SENSORGRAMS ---
% =========================================================================
fprintf('\n--- Step 6: Plotting results... ---\n');

% Create a time axis based on frame number and frame rate
time_axis = (0:numFrames-1) / frameRate;

figure('Name', 'Reconstructed Sensorgrams from Video', 'Position', [100, 100, 900, 600]);
plot(time_axis, ru_reconstructed, 'LineWidth', 1.5);

grid on;
box on;
title('Reconstructed Sensorgrams from Experimental Video');
xlabel('Time (s)');
ylabel('Response (RU)');
xlim([0, time_axis(end)]);

% Create a legend for each sensor line
legend_labels = arrayfun(@(j) sprintf('Line %d', j), 1:vidHeight, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northwest');

fprintf('Plotting complete. You can now analyze the `ru_reconstructed` data matrix.\n');

% Optional: Save the results to a .mat file for later use
% save('reconstructed_data.mat', 'time_axis', 'ru_reconstructed');
