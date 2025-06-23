% =========================================================================
% Main Script for RLS Parameter Estimation of Binding Kinetics
% Date: 2025-06-22
% =========================================================================
clear; close all; clc;

%% ------------------- SCRIPT CONFIGURATION ------------------------------
% Define the directory containing your Excel files
dir_path = 'C:\REPOS\PPBIO-Optical-Model\Adsorption\EXP_DATA\Data_HCL\data_excel';

% !!! CRUCIAL: Define the time (in seconds) when the concentration step occurs
% This is the time when you switch from buffer (C=0) to the analyte.
% This MUST be accurate for the estimation to work.
step_time = 0.0; % Example: Analyte introduced at 60 seconds

%% ------------------- DATA LOADING --------------------------------------
data_dict = load_experimental_data(dir_path);
file_names = fieldnames(data_dict);
n_experiments = length(file_names);
if n_experiments < 2
    error('This script requires at least two experimental files for initialization.');
end
fprintf('Loaded %d experimental files.\n', n_experiments);

%% ------------------- INITIALIZATION (NEW SECTION) ----------------------
% Get a data-driven initial guess for the parameters. This helps the RLS
% converge to a physically meaningful solution. We use the 2-point steady
% state solution from your original code.

fprintf('Calculating initial parameter guess...\n');
% Pass the filename string explicitly to the function
[R1, C1] = get_steady_state_point(data_dict.(file_names{1}), file_names{1}, step_time);
[R2, C2] = get_steady_state_point(data_dict.(file_names{2}), file_names{2}, step_time);
% Solve the 2x2 system for initial Ka and N_max
A = [R1, -C1; R2, -C2];
B = [-C1*R1; -C2*R2];
if det(A) == 0
    error('Cannot initialize parameters. The two steady-state points are not linearly independent.');
end
S = A \ B;
Ka_initial = abs(1/S(1));
N_max_initial = abs(S(2)); % Note: This is N_max, not N_max*theta_inf


% Convert physical guesses into an initial theta vector
gamma = -0;
beta = -gamma * abs(Ka_initial);
alpha = 1 + gamma / Ka_initial;
theta = [alpha; beta; gamma];

fprintf('Initial Guess -> Ka: %.2e, N_max: %.2e\n', Ka_initial, N_max_initial);

%% ------------------- RLS ESTIMATION ------------------------------------
% RLS Initial Parameters
P = eye(3) * 1e6;
lambda_ = 1.0;
theta_history = []; 

fprintf('Starting RLS estimation...\n');
for i = 1:n_experiments
    file_name = file_names{i};
    data = data_dict.(file_name);
    fprintf('Processing experiment: %s\n', file_name);

    % Pass the filename explicitly to the function
    [theta, P, history_for_exp] = run_rls_estimation(data, file_name, theta, P, lambda_, step_time);
    theta_history = [theta_history; history_for_exp];
end
fprintf('RLS estimation complete.\n');

%% ------------------- PARAMETER RECOVERY AND RESULTS --------------------
theta_final = theta; 
final_data = data_dict.(file_names{end});
time_s = final_data.Time / 1000;
dt_mean = mean(diff(time_s));

ka_final = -theta_final(3) / dt_mean;
kd_final = (1 - theta_final(1)) / dt_mean;
N_max_final = -theta_final(2) / theta_final(3);

fprintf('\n--- Final Estimated Physical Parameters ---\n');
fprintf('Association Rate (ka): %.4e\n', ka_final);
fprintf('Dissociation Rate (kd): %.4e\n', kd_final);
fprintf('Max. Capacity (N_max): %.4e\n', N_max_final);
fprintf('Affinity Constant (KA = ka/kd): %.4e\n', ka_final / kd_final);
fprintf('Dissociation Constant (KD = kd/ka): %.4e\n', kd_final / ka_final);

%% ------------------- PLOTTING AND VALIDATION ---------------------------
plot_parameter_convergence(theta_history, dt_mean);

% CRUCIAL FIX: Pass the file_names cell array to the plotting function
plot_model_validation(data_dict, file_names, ka_final, kd_final, N_max_final, step_time);


