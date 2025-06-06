% --- Main Script to Run the Iterative Experimental Design ---
clearvars; close all; clc;
fprintf('Starting Iterative Experimental Design Example (Focus on D-optimality & Condition Number)...\n');

% --- 1. Define Model Configuration ---
model_config.gridN_x = 10; 
model_config.gridN_y = 5; 
model_config.gridN_z = 3;
model_config.ads_layer = 1;
model_config.D_coeff = 6e-3;
model_config.ru_to_m = 1e-6;

% ** Using N_p = 6 (2 cells in adsorption region, y-direction heterogeneity) **
model_config.ads_x_range = [5,5]; 
model_config.ads_y_range = [2,3]; 
ads_param_shape = [model_config.ads_x_range(2)-model_config.ads_x_range(1)+1, ...
                   model_config.ads_y_range(2)-model_config.ads_y_range(1)+1]; % Should be [1,2]

% --- 2. Define Ground Truth Parameters for the Adsorption Region ---
rng(42); 
kon_base = 9.4e3; 
koff_base = 0.0078; 
num_cells_in_ads_region = prod(ads_param_shape);
smax_total_for_region = 1.0; 
smax_per_cell_base = smax_total_for_region / num_cells_in_ads_region;

true_kon_ads_region = kon_base * (1 + 0.1*randn(ads_param_shape)); 
true_koff_ads_region = koff_base * (1 + 0.1*randn(ads_param_shape));
true_smax_ads_region = smax_per_cell_base * (1 + 0.1*randn(ads_param_shape));

N_p_target = numel(true_kon_ads_region) + numel(true_koff_ads_region) + numel(true_smax_ads_region);
fprintf('Targeting N_p = %d parameters.\n', N_p_target);
if N_p_target == 0
    error('N_p_target is 0. Check adsorption region definition or parameter shapes.');
end

% --- 3. Create Candidate Experiments Pool (Revised for more orthogonality) ---
fprintf('Creating a revised candidate experiments pool...\n');
T_short = 300; T_med = 800; T_long = 2000; T_vlong_diss = 5000; 
C_vlow = 1e-7; C_low = 1e-6; C_mid = 5e-6; C_high = 2e-5; C_vhigh = 1e-4;
V_low = 2.0; V_mid = 8.3; V_high = 20.0; V_vhigh = 50.0;

candidate_experiments_pool = {};
base_s.c_diss = 0; 

% Set 1: koff Focus
s = base_s; s.pulse_times = [T_short]; s.pulse_concs = [C_high]; s.max_velocity = V_high; s.t_total = T_short + T_vlong_diss; candidate_experiments_pool{end+1}=s; % Exp C1
s = base_s; s.pulse_times = [T_med]; s.pulse_concs = [C_mid]; s.max_velocity = V_mid; s.t_total = T_med + T_long; candidate_experiments_pool{end+1}=s; % Exp C2

% Set 2: Smax / kon Focus (Saturation)
s = base_s; s.pulse_times = [T_long]; s.pulse_concs = [C_vhigh]; s.max_velocity = V_mid; s.t_total = T_long + T_short; candidate_experiments_pool{end+1}=s; % Exp C3
s = base_s; s.pulse_times = [T_med]; s.pulse_concs = [C_high]; s.max_velocity = V_high; s.t_total = T_med + T_med; candidate_experiments_pool{end+1}=s; % Exp C4

% Set 3: kon / Initial Rates / Concentration Dependence
s = base_s; s.pulse_times = [T_short]; s.pulse_concs = [C_low]; s.max_velocity = V_vhigh; s.t_total = T_short + T_med; candidate_experiments_pool{end+1}=s; % Exp C5
s = base_s; s.pulse_times = [T_short]; s.pulse_concs = [C_mid]; s.max_velocity = V_vhigh; s.t_total = T_short + T_med; candidate_experiments_pool{end+1}=s; % Exp C6
s = base_s; s.pulse_times = [T_short]; s.pulse_concs = [C_high]; s.max_velocity = V_vhigh; s.t_total = T_short + T_med; candidate_experiments_pool{end+1}=s; % Exp C7

% Set 4: Flow Rate Effect
s = base_s; s.pulse_times = [T_med]; s.pulse_concs = [C_mid]; s.max_velocity = V_low; s.t_total = T_med + T_long; candidate_experiments_pool{end+1}=s; % Exp C8
s = base_s; s.pulse_times = [T_med]; s.pulse_concs = [C_mid]; s.max_velocity = V_vhigh; s.t_total = T_med + T_long; candidate_experiments_pool{end+1}=s; % Exp C10 (Removed C9 as it was a duplicate of C2)

% Set 5: Different Pulse Durations
s = base_s; s.pulse_times = [T_short]; s.pulse_concs = [C_mid]; s.max_velocity = V_mid; s.t_total = T_short + T_long; candidate_experiments_pool{end+1}=s; % Exp C11
s = base_s; s.pulse_times = [T_long]; s.pulse_concs = [C_mid]; s.max_velocity = V_mid; s.t_total = T_long + T_long; candidate_experiments_pool{end+1}=s; % Exp C12

% Set 6: Multi-Step / Complex
T1_orig=800; T2_orig=1600; T3_orig=2400;
s = base_s; s.pulse_times = [T1_orig*0.5, T1_orig, T1_orig*1.5]; s.pulse_concs = [C_high, 0, C_mid]; s.max_velocity = V_mid; s.t_total = T1_orig*1.5 + T_med; candidate_experiments_pool{end+1}=s; % Exp C13
s = base_s; s.pulse_times = [T1_orig,T2_orig,T3_orig,T3_orig+T1_orig]; s.pulse_concs = [C_low,C_mid,C_high,0]; s.max_velocity = V_low; s.t_total = (T3_orig+T1_orig)+T_med; candidate_experiments_pool{end+1}=s; % Exp C14

% Set 7: Very Low Concentration
s = base_s; s.pulse_times = [T_long]; s.pulse_concs = [C_vlow]; s.max_velocity = V_mid; s.t_total = T_long + T_long; candidate_experiments_pool{end+1}=s; % Exp C15

fprintf('%d candidate experiments defined in the revised pool.\n', length(candidate_experiments_pool));

% --- 4. Call the Iterative Design Loop ---
max_exps = min(N_p_target + 4, length(candidate_experiments_pool)); 
if max_exps < N_p_target && N_p_target <= length(candidate_experiments_pool), max_exps = N_p_target; end
if max_exps == 0 && length(candidate_experiments_pool) > 0, max_exps = 1; end
if max_exps > length(candidate_experiments_pool), max_exps = length(candidate_experiments_pool); end

initial_selection = []; 

fprintf('Calling iterative_experimental_design_loop to select up to %d experiments...\n', max_exps);
[selected_settings, selected_indices, final_J, final_SVs, history] = ...
    iterative_experimental_design_loop(true_kon_ads_region, true_koff_ads_region, true_smax_ads_region, ...
                                     model_config, candidate_experiments_pool, N_p_target, max_exps, ...
                                     initial_selection, 1e-7, []); 

% --- DIAGNOSTIC BLOCK ---
if ~isempty(history)
    fprintf('\n--- DIAGNOSTIC: Fields in history(1) ---\n');
    disp(fieldnames(history(1)));
    if length(history) > 1 && length(history) == max_exps % Check last entry if loop completed
        fprintf('\n--- DIAGNOSTIC: Fields in history(end) ---\n');
        disp(fieldnames(history(end)));
    end
else
    fprintf('\n--- DIAGNOSTIC: History is empty! ---\n');
end
% --- END DIAGNOSTIC BLOCK ---

% --- 5. Analyze and Plot Results ---
if ~isempty(history)
    fprintf('\n--- OED Loop Results --- \n');
    for i=1:length(history)
        fprintf('Step %d: %d Exps, Indices: [%s], Rank: %d, MinSV_nz: %.2e, LogDet: %.2e, Cond: %.2e\n', ...
            i, history(i).num_experiments, num2str(history(i).selected_indices_so_far), ...
            history(i).rank, history(i).min_singular_value_nz, history(i).log_det_JtJ, history(i).condition_number);
    end

    figure('Name', 'OED History (D-Optimality Focus)');
    subplot(4,1,1); 
    plot([history.num_experiments], [history.rank], 'bo-'); 
    ylabel('Rank'); title(sprintf('Identifiability Improvement (N_p=%d)', N_p_target));
    ylim([0 N_p_target + 1]); grid on;

    subplot(4,1,2); 
    semilogy([history.num_experiments], [history.min_singular_value_nz], 'ro-'); 
    ylabel('Min Non-Zero SV (log)'); grid on;
    
    subplot(4,1,3);
    valid_logdet_idx = arrayfun(@(x) isfield(x, 'log_det_JtJ') && x.log_det_JtJ > -Inf, history);
    if any(valid_logdet_idx)
        plot([history(valid_logdet_idx).num_experiments], [history(valid_logdet_idx).log_det_JtJ], 'ms-');
    end
    ylabel('LogDet(JtJ)'); grid on;

    subplot(4,1,4);
    valid_cond_idx = arrayfun(@(x) isfield(x,'condition_number') && x.condition_number < Inf && x.condition_number > 0, history);
    if any(valid_cond_idx)
        semilogy([history(valid_cond_idx).num_experiments], [history(valid_cond_idx).condition_number], 'go-');
    end
    ylabel('Condition Number (log)'); xlabel('Number of Selected Experiments'); grid on;

    if ~isempty(final_SVs)
        final_hist_entry_idx = find([history.num_experiments] == length(selected_indices), 1, 'last');
        final_cond_num_to_display = Inf;
        if ~isempty(final_hist_entry_idx) && isfield(history(final_hist_entry_idx), 'condition_number')
            final_cond_num_to_display = history(final_hist_entry_idx).condition_number;
        end

        figure('Name', 'Final Singular Values (D-Optimality Focus)'); 
        semilogy(final_SVs, 'o-', 'LineWidth',1.5, 'MarkerSize',6); 
        title(sprintf('Final Singular Values for %d Selected Experiments (Rank %d, Cond %.2e)', ...
            length(selected_indices), history(end).rank, final_cond_num_to_display));
        ylabel('Singular Value Magnitude (log scale)'); xlabel('Singular Value Index');
        grid on;
        if history(end).rank < N_p_target && history(end).rank > 0 
            xline(history(end).rank + 0.5, 'r--', 'Label', sprintf('Rank Cutoff (%d)', history(end).rank));
        end
    end
    
    fprintf('\nSelected Experiment Settings:\n');
    if ~isempty(selected_settings)
        for i=1:length(selected_indices)
            current_setting = selected_settings{i};
            pulse_t_str = '[]'; if ~isempty(current_setting.pulse_times), pulse_t_str = num2str(current_setting.pulse_times(1)); end
            pulse_c_str = '[]'; if ~isempty(current_setting.pulse_concs), pulse_c_str = sprintf('%.1e',current_setting.pulse_concs(1)); end

            fprintf('Exp %d (Pool Index %d): Vel=%.1f, Conc1=%s, Npulses=%d, T_pulse1=%s, T_total=%.0f\n', ...
                i, selected_indices(i), current_setting.max_velocity, pulse_c_str, ...
                length(current_setting.pulse_times), pulse_t_str, current_setting.t_total);
        end
    end
else
    fprintf('The OED loop did not produce a history, check parameters or initial conditions.\n');
end



% --- Main Iterative Experimental Design Function (Modified for D-optimality & Condition Number) ---
function [selected_experiments_settings, selected_experiments_indices, final_Jacobian, final_singular_values, history] = ...
    iterative_experimental_design_loop(ground_truth_kon_ads_region, ground_truth_koff_ads_region, ground_truth_smax_ads_region, ...
                                     model_config, candidate_experiments_pool, N_p, max_experiments_to_select, ...
                                     initial_selected_indices, rank_tol_factor, min_sval_target_abs)

    fprintf('>>> EXECUTING OED LOOP VERSION WITH log_det_JtJ (Focus: D-Optimality) <<<\n');
    % --- Input Defaults & Initialization ---
    if ~exist('initial_selected_indices', 'var') || isempty(initial_selected_indices)
        initial_selected_indices = [];
    end
    if ~exist('rank_tol_factor', 'var') || isempty(rank_tol_factor)
        rank_tol_factor = 1e-6;
    end
    if ~exist('min_sval_target_abs', 'var') 
        min_sval_target_abs = []; 
    end

    selected_experiments_indices = reshape(initial_selected_indices, 1, []);
    num_candidate_experiments = length(candidate_experiments_pool);
    selected_mask = false(1, num_candidate_experiments);
    selected_mask(selected_experiments_indices) = true;

    all_experiment_jacobian_blocks = cell(1, num_candidate_experiments); 

    p_true_flat_region = [ground_truth_kon_ads_region(:); ground_truth_koff_ads_region(:); ground_truth_smax_ads_region(:)];
    if length(p_true_flat_region) ~= N_p
        error('N_p (%d) does not match the size of provided ground truth parameters for adsorption region (%d).', N_p, length(p_true_flat_region));
    end
    num_kon_params_region = numel(ground_truth_kon_ads_region);
    num_koff_params_region = numel(ground_truth_koff_ads_region);

    history = [];
    Combined_Jacobian_current = [];
    current_num_selected_experiments = length(selected_experiments_indices);

    % If there are initial experiments, build their combined Jacobian and assess
    if current_num_selected_experiments > 0
        temp_J_rows = cell(current_num_selected_experiments, 1);
        fprintf('Processing %d initial experiments...\n', current_num_selected_experiments);
        for i = 1:current_num_selected_experiments
            idx = selected_experiments_indices(i);
            fprintf('Calculating Jacobian for initial experiment (Pool Index: %d)...\n', idx);
            all_experiment_jacobian_blocks{idx} = compute_jacobian_for_single_experiment( ...
                candidate_experiments_pool{idx}, p_true_flat_region, ...
                num_kon_params_region, num_koff_params_region, model_config, N_p);
            temp_J_rows{i} = all_experiment_jacobian_blocks{idx};
        end
        Combined_Jacobian_current = vertcat(temp_J_rows{:});
        
        if ~isempty(Combined_Jacobian_current)
            svals_hist = svd(Combined_Jacobian_current, 'econ');
            abs_tol_rank_hist = determine_rank_tolerance(svals_hist, rank_tol_factor);
            rank_hist = sum(svals_hist > abs_tol_rank_hist);
            svals_nz_hist = svals_hist(svals_hist > abs_tol_rank_hist);
            min_sval_nz_hist = 0; if ~isempty(svals_nz_hist), min_sval_nz_hist = min(svals_nz_hist); end
            
            log_det_hist = -Inf; cond_hist = Inf;
            if rank_hist == N_p
                relevant_svals_hist = svals_hist(1:N_p);
                if ~any(relevant_svals_hist <= abs_tol_rank_hist) && all(relevant_svals_hist > 0) % Ensure positive before log
                    log_det_hist = 2 * sum(log(relevant_svals_hist));
                end
                if min_sval_nz_hist > 0 && ~isempty(svals_hist) && svals_hist(1)>0, cond_hist = svals_hist(1)/min_sval_nz_hist; end
            end
            
            history_entry.num_experiments = current_num_selected_experiments;
            history_entry.selected_indices_so_far = selected_experiments_indices;
            history_entry.rank = rank_hist;
            history_entry.min_singular_value_nz = min_sval_nz_hist;
            history_entry.log_det_JtJ = log_det_hist; % Assignment
            history_entry.condition_number = cond_hist;
            history_entry.all_singular_values = svals_hist;
            history = [history, history_entry];
            fprintf('Initial set: Rank=%d/%d, MinSV_nz=%.2e, LogDet=%.2e, Cond=%.2e\n', rank_hist, N_p, min_sval_nz_hist, log_det_hist, cond_hist);

            if rank_hist == N_p && (~isempty(min_sval_target_abs) && min_sval_nz_hist >= min_sval_target_abs)
                 fprintf('Initial set already meets rank and MinSV target. Stopping.\n');
                 selected_experiments_settings = candidate_experiments_pool(selected_experiments_indices);
                 final_Jacobian = Combined_Jacobian_current;
                 final_singular_values = svals_hist;
                 return;
            end
        end
    end

    % --- Main Iterative Loop ---
    for exp_count_to_select = (current_num_selected_experiments + 1) : max_experiments_to_select
        best_next_candidate_pool_idx = -1;
        
        if isempty(history)
            best_score_rank = -1;
            best_score_min_sval = -Inf;
            best_score_log_det = -Inf; % Initialize D-optimality score
            best_score_condition_number = Inf; 
        else
            current_hist = history(end);
            best_score_rank = current_hist.rank;
            best_score_min_sval = current_hist.min_singular_value_nz;
            if isfield(current_hist, 'log_det_JtJ') % Check if field exists
                best_score_log_det = current_hist.log_det_JtJ; 
            else
                best_score_log_det = -Inf; % Fallback
            end
            best_score_condition_number = current_hist.condition_number;

            if current_hist.rank < N_p 
                best_score_log_det = -Inf;
                best_score_condition_number = Inf;
            end
        end

        fprintf('\n--- Selecting Experiment %d of up to %d ---\n', exp_count_to_select, max_experiments_to_select);
        
        num_available_to_test = num_candidate_experiments - (exp_count_to_select - 1);
        if num_available_to_test <= 0 && exp_count_to_select <= max_experiments_to_select
             fprintf('No more unique candidate experiments available from the pool.\n');
             break;
        end
        
        candidate_eval_counter = 0;
        for cand_pool_idx = 1:num_candidate_experiments
            if selected_mask(cand_pool_idx)
                continue; 
            end
            candidate_eval_counter = candidate_eval_counter + 1;
            fprintf('Evaluating candidate %d (Pool Index: %d)... ', candidate_eval_counter, cand_pool_idx);

            if isempty(all_experiment_jacobian_blocks{cand_pool_idx})
                all_experiment_jacobian_blocks{cand_pool_idx} = compute_jacobian_for_single_experiment( ...
                    candidate_experiments_pool{cand_pool_idx}, p_true_flat_region, ...
                    num_kon_params_region, num_koff_params_region, model_config, N_p);
            end
            J_block_candidate = all_experiment_jacobian_blocks{cand_pool_idx};

            if isempty(J_block_candidate) || size(J_block_candidate,2) ~= N_p
                fprintf('Warning: Jacobian for candidate (Pool Index %d) is invalid. Skipping.\n', cand_pool_idx);
                continue;
            end

            if isempty(Combined_Jacobian_current)
                J_test = J_block_candidate;
            else
                J_test = [Combined_Jacobian_current; J_block_candidate];
            end
            
            svals_test = svd(J_test, 'econ');
            abs_tol_rank_test = determine_rank_tolerance(svals_test, rank_tol_factor);
            rank_test = sum(svals_test > abs_tol_rank_test);
            svals_nz_test = svals_test(svals_test > abs_tol_rank_test);
            min_sval_test_nz = 0; if ~isempty(svals_nz_test), min_sval_test_nz = min(svals_nz_test); end
            
            log_det_test = -Inf;
            if rank_test == N_p 
                relevant_svals_for_D = svals_test(1:N_p); 
                if ~any(relevant_svals_for_D <= abs_tol_rank_test) && all(relevant_svals_for_D > 0) % Check positive for log
                    log_det_test = 2 * sum(log(relevant_svals_for_D));
                end
            end
            
            cond_test = Inf; 
            if rank_test == N_p && min_sval_test_nz > 0 && ~isempty(svals_test) && svals_test(1) > 0
                cond_test = svals_test(1) / min_sval_test_nz;
            end
            fprintf('[R:%d,MinSVnz:%.2e,LogDet:%.2e,Cond:%.2e] ', rank_test, min_sval_test_nz, log_det_test, cond_test);

            improved = false;
            if rank_test > best_score_rank
                improved = true; 
            elseif rank_test == best_score_rank
                if rank_test == N_p 
                    if log_det_test > best_score_log_det 
                        improved = true;
                    elseif (abs(log_det_test - best_score_log_det) < 1e-1 || ... 
                           (best_score_log_det > -1e6 && abs(best_score_log_det)>1e-9 && abs(log_det_test - best_score_log_det) / abs(best_score_log_det) < 0.001))
                        if cond_test < best_score_condition_number 
                            improved = true;
                        end
                    end
                else 
                    if min_sval_test_nz > best_score_min_sval
                        improved = true;
                    end
                end
            end

            if improved
                best_score_rank = rank_test;
                best_score_min_sval = min_sval_test_nz;
                best_score_log_det = log_det_test;
                best_score_condition_number = cond_test;
                best_next_candidate_pool_idx = cand_pool_idx;
                fprintf('-> New Best Candidate!\n');
            else
                fprintf('\n');
            end
        end 

        if best_next_candidate_pool_idx == -1
            fprintf('No candidate experiment found that improves selection criteria further.\n');
            break; 
        end

        selected_experiments_indices = [selected_experiments_indices, best_next_candidate_pool_idx];
        selected_mask(best_next_candidate_pool_idx) = true;
        
        J_block_selected_new = all_experiment_jacobian_blocks{best_next_candidate_pool_idx};
        if isempty(Combined_Jacobian_current)
            Combined_Jacobian_current = J_block_selected_new;
        else
            Combined_Jacobian_current = [Combined_Jacobian_current; J_block_selected_new];
        end
        current_num_selected_experiments = length(selected_experiments_indices);
        fprintf('Selected experiment with Pool Index: %d. (Total selected: %d)\n', best_next_candidate_pool_idx, current_num_selected_experiments);

        svals_current_combined = svd(Combined_Jacobian_current, 'econ');
        abs_tol_rank_final = determine_rank_tolerance(svals_current_combined, rank_tol_factor);
        rank_final = sum(svals_current_combined > abs_tol_rank_final);
        svals_nz_final = svals_current_combined(svals_current_combined > abs_tol_rank_final);
        min_sval_nz_final = 0; if ~isempty(svals_nz_final), min_sval_nz_final = min(svals_nz_final); end
        
        log_det_final = -Inf; cond_num_final = Inf; 
        if rank_final == N_p
            relevant_svals_final = svals_current_combined(1:N_p);
            if ~any(relevant_svals_final <= abs_tol_rank_final) && all(relevant_svals_final > 0) % Check positive for log
                log_det_final = 2 * sum(log(relevant_svals_final));
            end
            if min_sval_nz_final > 0 && ~isempty(svals_current_combined) && svals_current_combined(1)>0, cond_num_final = svals_current_combined(1)/min_sval_nz_final; end
        end
        
        fprintf('  Combined Jacobian Size: %d x %d\n', size(Combined_Jacobian_current,1), size(Combined_Jacobian_current,2));
        fprintf('  Rank: %d / %d\n', rank_final, N_p);
        fprintf('  Smallest Non-Zero Singular Value: %.2e\n', min_sval_nz_final);
        fprintf('  LogDet(JtJ): %.2e\n', log_det_final);
        fprintf('  Condition Number (approx): %.2e\n', cond_num_final);

        history_entry.num_experiments = current_num_selected_experiments;
        history_entry.selected_indices_so_far = selected_experiments_indices;
        history_entry.rank = rank_final;
        history_entry.min_singular_value_nz = min_sval_nz_final;
        history_entry.log_det_JtJ = log_det_final; % Assignment
        history_entry.condition_number = cond_num_final;
        history_entry.all_singular_values = svals_current_combined;
        history = [history, history_entry];

        if rank_final == N_p
            if ~isempty(min_sval_target_abs) && min_sval_nz_final >= min_sval_target_abs
                fprintf('Target rank achieved AND minimum singular value target met. Stopping.\n');
                break;
            end
        end
        if current_num_selected_experiments >= max_experiments_to_select 
            fprintf('Reached maximum number of experiments to select. Stopping.\n');
            break;
        end
    end 

    selected_experiments_settings = candidate_experiments_pool(selected_experiments_indices);
    final_Jacobian = Combined_Jacobian_current;
    if ~isempty(final_Jacobian)
        final_singular_values = svd(final_Jacobian, 'econ');
    else
        final_singular_values = [];
    end

    fprintf('\n--- Iterative Design Loop Finished ---\n');
    fprintf('Selected a total of %d experiments.\n', length(selected_experiments_indices));
    fprintf('Selected experiment indices from pool: %s\n', mat2str(selected_experiments_indices));
    if ~isempty(history)
        final_hist_entry = history(end);
        fprintf('Final Rank: %d/%d, MinSV_nz: %.2e, LogDet: %.2e, Cond: %.2e\n', ...
            final_hist_entry.rank, N_p, final_hist_entry.min_singular_value_nz, ...
            final_hist_entry.log_det_JtJ, final_hist_entry.condition_number);
    else
        fprintf('No experiments were selected or processed.\n');
    end
end

% --- Helper function to compute Jacobian for a single experiment ---
function J_exp = compute_jacobian_for_single_experiment(experiment_setting, p_true_flat_region, ...
                                                    num_kon_region, num_koff_region, model_config, N_p_total)
    % Unpack model_config
    gridN_x = model_config.gridN_x; gridN_y = model_config.gridN_y; gridN_z = model_config.gridN_z;
    ads_x_range = model_config.ads_x_range; ads_y_range = model_config.ads_y_range; ads_layer = model_config.ads_layer;
    D_coeff = model_config.D_coeff; ru_to_m = model_config.ru_to_m;
    
    ads_region_shape_x = ads_x_range(2)-ads_x_range(1)+1;
    ads_region_shape_y = ads_y_range(2)-ads_y_range(1)+1;
    ads_param_shape = [ads_region_shape_x, ads_region_shape_y];

    kon_base_ads = reshape(p_true_flat_region(1:num_kon_region), ads_param_shape);
    koff_base_ads = reshape(p_true_flat_region(num_kon_region+1 : num_kon_region+num_koff_region), ads_param_shape);
    smax_base_ads = reshape(p_true_flat_region(num_kon_region+num_koff_region+1 : N_p_total), ads_param_shape); % Corrected end index

    [~, s_obs_base] = run_single_experiment_for_oed(gridN_x, gridN_y, gridN_z, ...
        kon_base_ads, koff_base_ads, smax_base_ads, ...
        ads_x_range, ads_y_range, ads_layer, experiment_setting, D_coeff, ru_to_m, model_config); % Pass model_config

    if isempty(s_obs_base)
        warning('Base sim for Jacobian (vel=%.1f, conc1=%.1e) returned empty s_obs.', experiment_setting.max_velocity, experiment_setting.pulse_concs(1));
        J_exp = zeros(0, N_p_total); 
        return;
    end
    num_time_points = length(s_obs_base);
    if num_time_points == 0 % Handle case where s_obs_base might be empty but not caught by isempty
        warning('Base sim for Jacobian (vel=%.1f, conc1=%.1e) resulted in 0 time points.', experiment_setting.max_velocity, experiment_setting.pulse_concs(1));
        J_exp = zeros(0, N_p_total);
        return;
    end
    J_exp = zeros(num_time_points, N_p_total);
    h_rel = 1e-5; 

    for k = 1:N_p_total
        p_pert_flat = p_true_flat_region;
        val_pk = p_true_flat_region(k);
        h_abs = h_rel * val_pk;
        if abs(val_pk) < 1e-9 || abs(h_abs) < 1e-12 % If param is zero or relative step is too small
             h_abs = h_rel; % Try a small absolute step based on h_rel as a magnitude
             if abs(val_pk) > 1 && h_abs > 1e-3 * abs(val_pk) % if val_pk is large, h_rel might be too large
                 h_abs = 1e-5 * abs(val_pk); % smaller relative step
             end
        end
        if abs(h_abs) < 1e-12, h_abs = 1e-8 * (1+abs(val_pk)); end % Final fallback, ensure it scales if val_pk is large
        if abs(h_abs) < 1e-12, h_abs = 1e-8; end


        p_pert_flat(k) = p_true_flat_region(k) + h_abs;

        kon_pert_ads = reshape(p_pert_flat(1:num_kon_region), ads_param_shape);
        koff_pert_ads = reshape(p_pert_flat(num_kon_region+1 : num_kon_region+num_koff_region), ads_param_shape);
        smax_pert_ads = reshape(p_pert_flat(num_kon_region+num_koff_region+1 : N_p_total), ads_param_shape); % Corrected end index

        [~, s_obs_pert] = run_single_experiment_for_oed(gridN_x, gridN_y, gridN_z, ...
            kon_pert_ads, koff_pert_ads, smax_pert_ads, ...
            ads_x_range, ads_y_range, ads_layer, experiment_setting, D_coeff, ru_to_m, model_config); % Pass model_config
        
        if length(s_obs_pert) ~= num_time_points
             warning('Jacobian calc: Length s_obs_pert (%d) ~= s_obs_base (%d) for param %d. Using NaNs. Vel=%.1f', length(s_obs_pert), num_time_points, k, experiment_setting.max_velocity);
             if num_time_points > 0, J_exp(:, k) = NaN; end % Fill with NaNs
             % Attempt to fix by truncating/padding if desperate, but it's risky
             % For now, NaNs will propagate and be handled.
             continue; 
        end
        
        if h_abs == 0 % Should be prevented by checks above
            J_exp(:,k) = 0;
        else
            J_exp(:, k) = (s_obs_pert - s_obs_base) / h_abs;
        end
    end
    
    if any(isinf(J_exp(:))) || any(isnan(J_exp(:)))
        warning('Jacobian (vel=%.1f) has NaNs/Infs. NaNs->0, Infs->big value.', experiment_setting.max_velocity);
        J_exp(isnan(J_exp)) = 0; 
        J_exp(isinf(J_exp) & J_exp > 0) = 1e12;
        J_exp(isinf(J_exp) & J_exp < 0) = -1e12;
    end
end

% --- Helper function to adapt from user's code for OED. Pass model_config for create_velocity_profile ---
function [t, s_obs] = run_single_experiment_for_oed(...
    nx, ny, nz, kon_ads, koff_ads, smax_ads,... 
    ads_x_range, ads_y_range, ads_layer, setting, D_coeff, ru_to_m, model_config) % Added model_config
    
    kon_grid = zeros(nx, ny, nz); koff_grid = zeros(nx, ny, nz); smax_grid = zeros(nx, ny, nz);
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_ads;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_ads;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_ads;
    
    % Ensure create_velocity_profile, simulate_3d_flow_model_with_pulses, compute_s_obs are in path
    % And that they are the versions provided by the user.
    velocity_profile = create_velocity_profile(nz, setting.max_velocity); % User's function
    t_breaks = [0, setting.pulse_times, setting.t_total];
    concentrations = [setting.pulse_concs, setting.c_diss];
    s0_grid = zeros(nx, ny, nz); % Initial surface concentration is zero for each simulation run here.
    
    % Assuming simulate_3d_flow_model_with_pulses is the user's version.
    % Its signature in user code: [t, c_s, s, K, Q] = simulate_3d_flow_model_with_pulses(...)
    [t, ~, s_grid_all_time, ~, ~] = simulate_3d_flow_model_with_pulses(...
        nx, ny, nz, kon_grid, koff_grid, smax_grid, ...
        velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid);
    
    s_obs = compute_s_obs(s_grid_all_time, ads_x_range, ads_y_range, ads_layer); % User's function
end

% --- Helper for rank tolerance ---
function tol = determine_rank_tolerance(svals, factor)
    if isempty(svals) || svals(1) == 0 % Handle case of zero matrix or all zero svals
        tol = 1e-9; 
    else
        tol = factor * svals(1); 
        if tol < 1e-12 
            tol = 1e-12;
        end
    end
end

function [kon_grid, koff_grid, smax_grid] = create_heterogeneous_grids_from_ads(...
    nx, ny, nz, ads_x_range, ads_y_range, ads_layer, kon_ads, koff_ads, smax_ads)
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_ads;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_ads;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_ads;
end

function s_obs = compute_s_obs(s_grid, ads_x_range, ads_y_range, ads_layer)
    ads_cells = s_grid(:, ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer);
    s_obs = squeeze(sum(ads_cells, [2,3,4]));
end

function [kon_grid, koff_grid, smax_grid] = create_homogeneous_grids(params, nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    kon = params(1);
    koff = params(2);
    smax_total = params(3);
    smax_per_cell = smax_total / ((ads_x_range(2)-ads_x_range(1)+1)*(ads_y_range(2)-ads_y_range(1)+1));
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_per_cell;
end

function velocity_profile = create_velocity_profile(nz, max_velocity)
    z_indices = 0:(nz-1);
    h = nz-1;
    velocity_profile = 4 * max_velocity * (z_indices/h) .* (1 - z_indices/h);
    velocity_profile = reshape(velocity_profile, [1,1,nz]);
end

function [kon_grid, koff_grid, smax_grid] = create_ground_truth_heterogeneity(nx, ny, nz, ads_x_range, ads_y_range, ads_layer)
    kon_base = 9.4e3;
    koff_base = 0.0078;
    smax_total = 1.0;
    num_ads_cells = (ads_x_range(2)-ads_x_range(1)+1) * (ads_y_range(2)-ads_y_range(1)+1);
    
    kon_grid = zeros(nx, ny, nz);
    koff_grid = zeros(nx, ny, nz);
    smax_grid = zeros(nx, ny, nz);
    
    rng(42);
    kon_vals = kon_base * (1 + 0.05*randn(ads_x_range(2)-ads_x_range(1)+1, ads_y_range(2)-ads_y_range(1)+1));
    koff_vals = koff_base * (1 + 0.05*randn(size(kon_vals)));
    smax_vals = (smax_total/num_ads_cells) * (1 + 0.05*randn(size(kon_vals)));
    
    kon_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = kon_vals;
    koff_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = koff_vals;
    smax_grid(ads_x_range(1):ads_x_range(2), ads_y_range(1):ads_y_range(2), ads_layer) = smax_vals;
end

function [t, c_s, s, K, Q] = simulate_3d_flow_model_with_pulses(...
    nx, ny, nz, kon_grid, koff_grid, smax_grid, velocity_profile, t_breaks, concentrations, D_coeff, ru_to_m, s0_grid)
    
    % Initialize state variables
    num_cells = nx * ny * nz;
    c_s0 = zeros(nx, ny, nz);
    c_s0(1, :, :) = concentrations(1); % Initial concentration
    s0 = s0_grid;
    Q0 = zeros(nx, ny, nz);
    R0 = zeros(nx, ny, nz);
    y0 = [c_s0(:); s0(:); Q0(:); R0(:)];
    
    % Setup ODE options
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    
    % Preallocate results
    t_all = [];
    y_all = [];
    
    % Process each time segment
    num_segments = length(t_breaks) - 1;
    for seg = 1:num_segments
        t_start = t_breaks(seg);
        t_end = t_breaks(seg+1);
        c0_seg = concentrations(seg);
        
        % Determine time points for segment
        num_points = max(10, ceil(200 * (t_end - t_start) / (t_breaks(end) - t_breaks(1))));
        tspan = linspace(t_start, t_end, num_points);
        
        % Run simulation for segment
        [t_seg, y_seg] = ode15s(@(t,y) ode_system(t, y, nx, ny, nz, velocity_profile, ...
            kon_grid, koff_grid, smax_grid, c0_seg, D_coeff, ru_to_m), tspan, y0, options);
        
        % Handle first segment specially
        if seg == 1
            t_all = t_seg;
            y_all = y_seg;
        else
            % Append results (skip first point to avoid duplicate)
            t_all = [t_all; t_seg(2:end)];
            y_all = [y_all; y_seg(2:end, :)];
        end
        
        % Update initial condition for next segment
        if seg < num_segments
            y0 = y_seg(end, :)';
            c_s_end = reshape(y0(1:num_cells), [nx, ny, nz]);
            s_end = reshape(y0(num_cells+1:2*num_cells), [nx, ny, nz]);
            c_s_end(1, :, :) = concentrations(seg+1);
            y0 = [c_s_end(:); s_end(:); y0(2*num_cells+1:end)]; % Preserve Q/R
        end
    end
    
    % Extract variables
    c_s = reshape(y_all(:, 1:num_cells), [length(t_all), nx, ny, nz]);
    s = reshape(y_all(:, num_cells+1:2*num_cells), [length(t_all), nx, ny, nz]);
    Q = reshape(y_all(:, 2*num_cells+1:3*num_cells), [length(t_all), nx, ny, nz]);
    R = reshape(y_all(:, 3*num_cells+1:end), [length(t_all), nx, ny, nz]);
    
    % Compute kernel
    K = exp(-Q) .* R;
    t = t_all;
end

function dydt = ode_system(t, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0, D_coeff, ru_to_m)
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);
    Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);
    R = reshape(y(3*num_cells + 1:4*num_cells), [nx, ny, nz]);
    dcsdt = zeros(nx, ny, nz);
    dsdt = zeros(nx, ny, nz);

    % Diffusion terms
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = (c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:));
    
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2);
    d2c_dz2(:,:,1) = c_s(:,:,2) - 2*c_s(:,:,1) + c_s(:,:,1);
    d2c_dz2(:,:,end) = c_s(:,:,end-1) - 2*c_s(:,:,end) + c_s(:,:,end-1);
    
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2);
    
    % Advection
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)));
    
    % Adsorption kinetics
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c0;
    dcsdt(1,:,:) = 0;
    
    % Compute dQ/dt and dR/dt
    dQdt = kon_grid .* c_s + koff_grid;
    dRdt = c_s .* exp(Q);

    % Combine all derivatives
    dydt = [dcsdt(:); dsdt(:); dQdt(:); dRdt(:)];
end