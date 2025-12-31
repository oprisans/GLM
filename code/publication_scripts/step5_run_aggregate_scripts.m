% orchestrator_run_minimal.m
% Minimal pipeline orchestrator: set grid, export to base, verify steps, run.
% Edit the grid below as needed.

%% 0) grid / seeds (edit here)
seedList     = 1:30;
Nvals        = [10 20 30];
sparsityList = 0.05:0.05:0.50;
snr_vals     = [0 5 10];
Tval         = 12000;

% expose to base workspace for scripts that read base variables
assignin('base','seedList', seedList);
assignin('base','Nvals', Nvals);
assignin('base','sparsityList', sparsityList);
assignin('base','snr_vals', snr_vals);
assignin('base','Tval', Tval);

%% 1) steps list (exact filenames)
steps = {
 'step1_driver_simulate_var_calcium.m', ...
 'step2_aa1_lassoglm_1y_ALLx_zScore.m', ...
 'step3_compare_groundTruth_vs_aa1_lassoglm_1x_ALLy_zScore.m', ...
 'step4_aggregate_and_plot_zScore.m' };

%% 2) resolve and run steps
for i = 1:numel(steps)
    sname = steps{i};
    fprintf('\n--- STEP %d: %s ---\n', i, sname);

    % try which() first, then fallback to cwd
    sp = which(sname);
    if isempty(sp)
        cand = fullfile(pwd, sname);
        if exist(cand,'file')==2
            sp = cand;
        end
    end

    if isempty(sp)
        fprintf('MISSING: %s (not on path and not in cwd) - skipping\n', sname);
        continue;
    end

    % Special-case: prefer calling aggregate_and_plot_zScore as a function
    if strcmpi(sname, 'step4_aggregate_and_plot_zScore.m')
        % if function exists on path call it directly (safe, gives correct args)
        if exist('aggregate_and_plot_zScore','file')==2 || exist('aggregate_and_plot_zScore','builtin')==5
            try
                fprintf('Calling function: aggregate_and_plot_zScore\n');
                aggregate_and_plot_zScore({'auprc_pos','F1_at_match'}, Nvals, true);
                fprintf('STEP %d completed: %s (function)\n', i, sname);
            catch ME
                warning('STEP %d FAILED when calling function aggregate_and_plot_zScore: %s', i, ME.message);
            end
        else
            % fallback: run the script file
            try
                run(sp);
                fprintf('STEP %d completed: %s (ran script)\n', i, sname);
            catch ME
                warning('STEP %d FAILED running script %s: %s', i, sp, ME.message);
            end
        end
    else
        % normal step: run the script (will execute in caller workspace)
        try
            run(sp);
            fprintf('STEP %d completed: %s\n', i, sname);
        catch ME
            warning('STEP %d FAILED running script %s: %s', i, sp, ME.message);
        end
    end
end

fprintf('\nMinimal orchestration finished.\n');