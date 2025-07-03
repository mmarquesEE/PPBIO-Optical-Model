clearvars; close all; clc;
rng(42); % For reproducibility

%% PART 1: SETUP THE PHYSICAL MODEL
fprintf('Setting up physical model parameters...\n');
% --- Define Model and Experiment Parameters ---
gridN_x = 22; gridN_y = 5; gridN_z = 3; ads_layer = 1;
ads_x_range = [5,15]; ads_y_range = [1,5];
ads_y_dim = ads_y_range(2) - ads_y_range(1) + 1;
D_coeff = 6e-5; ru_to_m = 1e-10; base_max_velocity = 8.3;
grid_size_x = 11.0; grid_size_z = 0.3;
dx = grid_size_x / gridN_x; dz = grid_size_z / gridN_z;
n_exp = 1; T1 = 2200; T2 = 2*T1; T3 = 3*T1; t_total = 4*T1;
c_diss = 0; c1 = 3.3e-4; c2 = 0.21e-4;

exp_settings = generate_experiments(n_exp, base_max_velocity, T1, T2, T3, t_total, c_diss, c1, c2);
model_config.gridN_x = gridN_x; model_config.gridN_y = gridN_y; model_config.gridN_z = gridN_z;
model_config.dx = dx; model_config.dz = dz;
model_config.ads_x_range = ads_x_range; model_config.ads_y_range = ads_y_range;
model_config.ads_layer = ads_layer;
model_config.D_coeff = D_coeff; model_config.ru_to_m = ru_to_m;

% We need the dimensions for the network output layer, so we run the
% simulation once to get the size of the sensorgram.
p_temp = generate_random_params(ads_y_dim, struct('kon',[3,4],'koff',[-3,-2],'smax',[2,3]), 1);
[t_eval, s_temp] = run_single_experiment_1D_model(p_temp, exp_settings(1), model_config);
fprintf('Model setup complete.\n\n');

%% PART 2: DEFINE THE SURROGATE NEURAL NETWORK
fprintf('Defining surrogate network architecture...\n');
num_params = 3 * ads_y_dim;
num_outputs = numel(s_temp);

% A deeper network is better for learning complex physics
layers = [
    featureInputLayer(num_params, 'Name', 'params')
    fullyConnectedLayer(256, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(512, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(256, 'Name', 'fc3')
    reluLayer('Name', 'relu3')
    fullyConnectedLayer(num_outputs, 'Name', 'output')
    ];
dlnet = dlnetwork(layers);
fprintf('Network defined.\n\n');

%% PART 3: TRAIN THE SURROGATE MODEL
% --- Training Setup ---
epochs = 500;      % RECOMMENDATION: Increase to 500-2000 for high accuracy
batchSize = 8;      % Number of simulations to run in parallel per epoch
learningRate = 0.001;

trailingAvg = [];
trailingAvgSq = [];

% Define reasonable bounds for generating random kinetic parameters
bounds.kon  = [3, 6];
bounds.koff = [-4, -1];
bounds.smax = [1, 3.5];

fprintf('Starting surrogate model training...\n');
fprintf('This will call the full physics model %d times.\n', epochs * batchSize);
tic;

lossHistory = [];
for epoch = 1:epochs
    p_batch = generate_random_params(ads_y_dim, bounds, batchSize);

    s_physics_batch = zeros(num_outputs, batchSize);
    % Use parfor for parallel execution to speed up data generation
    parfor i = 1:batchSize
        [~, s_physics] = run_single_experiment_1D_model(p_batch(:, i), exp_settings(1), model_config);
        s_physics_batch(:, i) = s_physics(:);
    end
    
    % =============================== FIX IS HERE ===============================
    % Manually filter out any columns corresponding to failed simulations (NaNs)
    valid_cols = ~any(isnan(s_physics_batch), 1);
    if ~any(valid_cols)
        fprintf('Epoch %d/%d, Skipped - all physics simulations in batch failed.\n', epoch, epochs);
        continue;
    end
    
    p_batch_clean = p_batch(:, valid_cols);
    s_physics_batch_clean = s_physics_batch(:, valid_cols);
    
    % Convert the clean batch to dlarray for training
    dl_p_batch = dlarray(p_batch_clean, 'CB');
    dl_s_physics_batch = dlarray(s_physics_batch_clean, 'CB');
    % ===========================================================================

    if canUseGPU
        dl_p_batch = gpuArray(dl_p_batch);
        dl_s_physics_batch = gpuArray(dl_s_physics_batch);
    end

    [loss, gradients] = dlfeval(@modelLoss, dlnet, dl_p_batch, dl_s_physics_batch);

    [dlnet, trailingAvg, trailingAvgSq] = adamupdate(dlnet, gradients, trailingAvg, trailingAvgSq, epoch, learningRate);

    lossHistory = [lossHistory, extractdata(loss)];
    if mod(epoch, 20) == 0 || epoch == 1
      fprintf('Epoch %d/%d, Loss: %.4e\n', epoch, epochs, loss);
    end
end
trainingTime = toc;
fprintf('Surrogate model training complete. Time taken: %.2f seconds.\n\n', trainingTime);

%% PART 4: VALIDATE THE TRAINED SURROGATE MODEL
fprintf('--- Validating Surrogate Model Performance ---\n\n');

% --- 1. Generate a new, unseen test parameter set ---
p_test = generate_random_params(ads_y_dim, bounds, 1);
fprintf('Generated new test parameters.\n');

% --- 2. Time the original, slow physics simulation ---
fprintf('Running original physics model for timing and accuracy comparison...\n');
tic;
[~, s_physics_test] = run_single_experiment_1D_model(p_test, exp_settings(1), model_config);
time_physics = toc;
fprintf('Original Physics Model execution time: %.4f seconds.\n', time_physics);

% --- 3. Time the trained, fast surrogate network ---
fprintf('Running trained surrogate model for timing comparison...\n');
dl_p_test = dlarray(p_test, 'CB');
if canUseGPU
    dlnet = dlupdate(@gpuArray, dlnet);
    dl_p_test = gpuArray(dl_p_test);
end
tic;
s_surrogate_pred_flat = predict(dlnet, dl_p_test);
time_surrogate = toc;
fprintf('Trained Surrogate Model execution time: %.4f seconds.\n', time_surrogate);

% --- 4. Compare Performance ---
speedup_factor = time_physics / time_surrogate;
fprintf('\n>>> SPEEDUP: The trained surrogate is approximately %.0fx faster. <<<\n\n', speedup_factor);

% --- 5. Compare Accuracy ---
s_surrogate_test = reshape(extractdata(gather(s_surrogate_pred_flat)), size(s_temp));
accuracy_rmse = sqrt(mean((s_physics_test - s_surrogate_test).^2, 'all'));
fprintf('>>> ACCURACY: Root Mean Squared Error (RMSE) between surrogate and physics model is %.4e.\n', accuracy_rmse);

% --- 6. Visualize the Comparison ---
figure('Position', [100, 100, 1200, 600]);
sgtitle('Surrogate Model Validation on Unseen Test Data');

subplot(1, 2, 1);
plot(lossHistory);
title('Surrogate Training Loss');
xlabel('Epoch');
ylabel('MSE Loss');
grid on;
set(gca, 'YScale', 'log');

subplot(1, 2, 2);
plot(t_eval, s_physics_test(:, 1), 'b-', 'LineWidth', 3, 'DisplayName', 'Original Physics Model');
hold on;
plot(t_eval, s_surrogate_test(:, 1), 'r--', 'LineWidth', 2, 'DisplayName', 'Trained Surrogate');
hold off;
title('Accuracy Comparison (First Sensorgram Line)');
xlabel('Time (s)');
ylabel('Response (RU)');
legend;
grid on;
xlim([0, t_eval(end)]);

% At the end, you can save your trained model for future use
% save('fast_forward_model.mat', 'dlnet');

%% HELPER FUNCTIONS
function [loss, gradients] = modelLoss(dlnet, p_batch, s_physics_batch)
    s_pred_batch = predict(dlnet, p_batch);
    % =============================== FIX IS HERE ===============================
    % The 'MissingData' option is removed to support older MATLAB versions.
    % NaN filtering is now handled in the main training loop.
    loss = mse(s_pred_batch, s_physics_batch);
    % ===========================================================================
    gradients = dlgradient(loss, dlnet.Learnables);
end

function p_batch = generate_random_params(num_lines, bounds, batchSize)
    num_params_per_line = 3;
    p_batch = zeros(num_lines * num_params_per_line, batchSize);
    for i = 1:batchSize
        log_kon  = bounds.kon(1)  + (bounds.kon(2)-bounds.kon(1))   * rand(num_lines, 1);
        log_koff = bounds.koff(1) + (bounds.koff(2)-bounds.koff(1)) * rand(num_lines, 1);
        log_smax = bounds.smax(1) + (bounds.smax(2)-bounds.smax(1)) * rand(num_lines, 1);
        p_batch(:, i) = [10.^log_kon; 10.^log_koff; 10.^log_smax];
    end
end
