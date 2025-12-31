% compare_groundtruth_vs_glm_full.m
% Complete comparison between ground-truth A (sheet 'A') and GLM result b1_CV1
% Computes: AUPRC, AUROC (pos/neg), Pearson r (pos/neg/overall), MAE, RMSE,
% sign accuracy, density, SHD, signed-SHD, weight correlations, and bootstrap CIs (optional).
%
% USAGE:
%  - Edit the filename construction block below to match your file naming.
%  - Set nboot (0 to skip bootstrap).
%  - Run the script.

clear all; close all; clc;
fullpath = mfilename('fullpath');
if isempty(fullpath)
    fullpath = matlab.desktop.editor.getActiveFilename;
end
[scriptFolder, scriptName, ~] = fileparts(fullpath);
parentFolder = fileparts(scriptFolder);

dataPath   = fullfile(parentFolder, 'data');
resultsPath= fullfile(parentFolder, 'results');

%% ---------- user config (match your filenames) ----------

folderA='groundTruth';
pathIn_A = fullfile(dataPath,folderA);

GLM_method ='lassoglm';
dataType='zScore';
folderB=extractAfter(scriptName,'compare_groundTruth_vs_');
pathIn_B = fullfile(dataPath,dataType,folderB);

pathOut = fullfile(resultsPath,dataType,folderB);
if ~exist(pathOut , 'dir')
    mkdir(pathOut );
end

% parameters used to build filenames (must match how files were saved)
Nval = 10;
Tval = 12000;
dtval = 0.05;
sparsityval = 0.05;
snr_dbval = 0;
seedval = 42;
wval = 0.25;
sigma_rval = 0.2;
tauval = 1.2;
rho_targetval = 0.9;
smooth_sigmaval = 0.2;

% how many bootstrap samples for CI (0 to disable)
nboot = 200;

%% ---------- build filename (keeps same exact formatting you used) ----------
  idx=1;
 %for Nval=10:10:30;
  %   for Tval=2400:1200:12000
         for sparsityval=0.05:0.05:0.3
             for snr_dbval=0:5:10
                  fileName=['N_',num2str(Nval),'_T_',num2str(Tval),...
                                        '_dt_',num2str(dtval),'_sparsity_',num2str(sparsityval),...
                                        '_snr_',num2str(snr_dbval),'_sigma_',num2str(sigma_rval),...
                                        '_tau_',num2str(tauval),'_w_',num2str(wval),'_rhoTarget_',...
                                        num2str(rho_targetval),'_smoothSigma_',num2str(smooth_sigmaval),...
                                        '_seed_',num2str(seedval),'.xlsx'];
               A_true = readmatrix(fullfile(pathIn_A,[folderA,'_',fileName]),'Sheet','A');
               b1_CV1 = readmatrix(fullfile(pathIn_B,[folderB,'_',fileName]),'Sheet','b1_CV1');
               b1_CV2 = readmatrix(fullfile(pathIn_B,[folderB,'_',fileName]),'Sheet','b1_CV2');
               b1_AIC = readmatrix(fullfile(pathIn_B,[folderB,'_',fileName]),'Sheet','b1_AIC');
               b1_BIC = readmatrix(fullfile(pathIn_B,[folderB,'_',fileName]),'Sheet','b1_BIC');
               b1_Custom = readmatrix(fullfile(pathIn_B,[folderB,'_',fileName]),'Sheet','b1_Custom');

                 %% ---------- compute metrics ----------
                results = compare_connectivity_matrices(A_true, b1_CV1, nboot);
                res_out_b1_CV1(idx,:)=[Nval, Tval, dtval, sparsityval, snr_dbval, seedval, wval, sigma_rval, tauval, ...
                    rho_targetval, smooth_sigmaval, results.auprc_pos, results.auprc_neg, ...
                    results.auroc_pos, results.auroc_neg, results.auprc_pos_CI, results.auprc_neg_CI, ...
                    results.r_pos, results.r_neg, results.r_overall, results.r_pred_on_predicted_edges, ...
                    results.mae_pos, results.mae_neg, results.mae_overall, results.rmse_pos, ...
                    results.rmse_neg, results.rmse_overall, results.sign_acc_pos, results.sign_acc_neg,...
                    results.sign_acc_overall, results.true_density, results.pred_density, results.SHD,...
                    results.SHD_signed, results.precision_at_match, results.recall_at_match, results.F1_at_match];

                results = compare_connectivity_matrices(A_true, b1_CV2, nboot);
                res_out_b1_CV2(idx,:)=[Nval, Tval, dtval, sparsityval, snr_dbval, seedval, wval, sigma_rval, tauval, ...
                    rho_targetval, smooth_sigmaval, results.auprc_pos, results.auprc_neg, ...
                    results.auroc_pos, results.auroc_neg, results.auprc_pos_CI, results.auprc_neg_CI, ...
                    results.r_pos, results.r_neg, results.r_overall, results.r_pred_on_predicted_edges, ...
                    results.mae_pos, results.mae_neg, results.mae_overall, results.rmse_pos, ...
                    results.rmse_neg, results.rmse_overall, results.sign_acc_pos, results.sign_acc_neg,...
                    results.sign_acc_overall, results.true_density, results.pred_density, results.SHD,...
                    results.SHD_signed, results.precision_at_match, results.recall_at_match, results.F1_at_match];

                results = compare_connectivity_matrices(A_true, b1_AIC, nboot);
                res_out_b1_AIC(idx,:)=[Nval, Tval, dtval, sparsityval, snr_dbval, seedval, wval, sigma_rval, tauval, ...
                    rho_targetval, smooth_sigmaval, results.auprc_pos, results.auprc_neg, ...
                    results.auroc_pos, results.auroc_neg, results.auprc_pos_CI, results.auprc_neg_CI, ...
                    results.r_pos, results.r_neg, results.r_overall, results.r_pred_on_predicted_edges, ...
                    results.mae_pos, results.mae_neg, results.mae_overall, results.rmse_pos, ...
                    results.rmse_neg, results.rmse_overall, results.sign_acc_pos, results.sign_acc_neg,...
                    results.sign_acc_overall, results.true_density, results.pred_density, results.SHD,...
                    results.SHD_signed, results.precision_at_match, results.recall_at_match, results.F1_at_match];

                results = compare_connectivity_matrices(A_true, b1_BIC, nboot);
                res_out_b1_BIC(idx,:)=[Nval, Tval, dtval, sparsityval, snr_dbval, seedval, wval, sigma_rval, tauval, ...
                    rho_targetval, smooth_sigmaval, results.auprc_pos, results.auprc_neg, ...
                    results.auroc_pos, results.auroc_neg, results.auprc_pos_CI, results.auprc_neg_CI, ...
                    results.r_pos, results.r_neg, results.r_overall, results.r_pred_on_predicted_edges, ...
                    results.mae_pos, results.mae_neg, results.mae_overall, results.rmse_pos, ...
                    results.rmse_neg, results.rmse_overall, results.sign_acc_pos, results.sign_acc_neg,...
                    results.sign_acc_overall, results.true_density, results.pred_density, results.SHD,...
                    results.SHD_signed, results.precision_at_match, results.recall_at_match, results.F1_at_match];
             
                results = compare_connectivity_matrices(A_true, b1_Custom, nboot);
                res_out_b1_Custom(idx,:)=[Nval, Tval, dtval, sparsityval, snr_dbval, seedval, wval, sigma_rval, tauval, ...
                    rho_targetval, smooth_sigmaval, results.auprc_pos, results.auprc_neg, ...
                    results.auroc_pos, results.auroc_neg, results.auprc_pos_CI, results.auprc_neg_CI, ...
                    results.r_pos, results.r_neg, results.r_overall, results.r_pred_on_predicted_edges, ...
                    results.mae_pos, results.mae_neg, results.mae_overall, results.rmse_pos, ...
                    results.rmse_neg, results.rmse_overall, results.sign_acc_pos, results.sign_acc_neg,...
                    results.sign_acc_overall, results.true_density, results.pred_density, results.SHD,...
                    results.SHD_signed, results.precision_at_match, results.recall_at_match, results.F1_at_match];
                idx=idx+1;

             end
         end
    % end
 %end

% Define variable names
variableNames = {'N', 'T', 'dt', 'sparsity', 'snr_db', 'seed', 'w', 'sigma_r', 'tau', 'rho_target', ...
    'smooth_sigma', 'auprc_pos', 'auprc_neg','auroc_pos', 'auroc_neg', 'auprc_pos_CI_low', ...
    'auprc_pos_CI_high','auprc_neg_CI_low', 'auprc_neg_CI_high','r_pos', 'r_neg', 'r_overall', 'r_pred_on_predicted_edges', 'mae_pos', ...
    'mae_neg', 'mae_overall', 'rmse_pos', 'rmse_neg', 'rmse_overall','sign_acc_pos', ...
    'sign_acc_neg', 'sign_acc_overall', 'true_density', 'pred_density', 'SHD', 'SHD_signed', ...
    'precision_at_match', 'recall_at_match', 'F1_at_match'};

    % Convert array to table with specified variable names
    T_b1_CV1 = array2table(res_out_b1_CV1, 'VariableNames', variableNames);
    % Write the table to the Excel file
    writetable(T_b1_CV1, [pathOut,filesep,[folderB,'_comp_groundTruth.xlsx']],'Sheet','b1_CV1');
    
    T_b1_CV2= array2table(res_out_b1_CV2, 'VariableNames', variableNames);
    % Write the table to the Excel file
    writetable(T_b1_CV2, [pathOut,filesep,[folderB,'_comp_groundTruth.xlsx']],'Sheet','b1_CV2');
    
    T_b1_AIC = array2table(res_out_b1_AIC, 'VariableNames', variableNames);
    % Write the table to the Excel file
    writetable(T_b1_AIC, [pathOut,filesep,[folderB,'_comp_groundTruth.xlsx']],'Sheet','b1_AIC');
    
    T_b1_BIC = array2table(res_out_b1_BIC, 'VariableNames', variableNames);
    % Write the table to the Excel file
    writetable(T_b1_BIC, [pathOut,filesep,[folderB,'_comp_groundTruth.xlsx']],'Sheet','b1_BIC');
    
    T_b1_Custom = array2table(res_out_b1_Custom, 'VariableNames', variableNames);
    % Write the table to the Excel file
    writetable(T_b1_Custom, [pathOut,filesep,[folderB,'_comp_groundTruth.xlsx']],'Sheet','b1_Custom');

save(fullfile(pathOut,[folderB,'_comp_groundTruth_b1_CV1']),"T_b1_CV1");
save(fullfile(pathOut,[folderB,'_comp_groundTruth_b1_CV2']),"T_b1_CV2");
save(fullfile(pathOut,[folderB,'_comp_groundTruth_b1_AIC']),"T_b1_AIC");
save(fullfile(pathOut,[folderB,'_comp_groundTruth_b1_BIC']),"T_b1_BIC");
save(fullfile(pathOut,[folderB,'_comp_groundTruth_b1_Custom']),"T_b1_Custom");

 % Nval=10:10:30; %col1
 % Tval=2400:1200:12000; %col2
 % sparsityval=0.05:0.05:0.3; %col4
 % snr_dbval=0:5:10; %col5

% for constant1 = 10:10:30;
%     constant4 = 0.05;
%     constant5 = 0;
% 
%     % --- Apply logical indexing to filter the data ---
%     % Create a logical array (a mask) where the conditions are met
%     % The & operator performs a logical AND operation
%     conditions = (res_out_b1_CV1(:, 1) == constant1) & (res_out_b1_CV1(:, 4) == constant4) & (res_out_b1_CV1(:, 5) == constant5);
%     % Extract the filtered a1 and y1 values
%     % A(conditions, 1) selects all rows where 'conditions' is true, and the first column (a1)
%     % A(conditions, end) selects all rows where 'conditions' is true, and the last column (y1)
%     a1_filtered = res_out_b1_CV1(conditions, 2);
%     y1_filtered = res_out_b1_CV1(conditions, 12); 
% 
%     % --- Plot the filtered data ---
%     %figure; % Open a new figure window
%     plot(a1_filtered, y1_filtered, 'o-'); % Plot y1 vs a1 with markers and a line
%     hold on
% 
%     % title(sprintf('y1 vs a1 for a2 = %d and a3 = %d', constant1, constant2));
%     % grid on;
% end
% xlabel('T');
% ylabel('auprc_pos');
% legend('10','20','30')


% %% ---------- display / table ----------
% disp('--- Summary results struct ---');
% disp(results.summary_table);
% 
% % print a short readable summary
% fprintf('\nQuick summary:\n');
% fprintf('AUPRC pos: %.4f   AUROC pos: %.4f\n', results.auprc_pos, results.auroc_pos);
% fprintf('AUPRC neg: %.4f   AUROC neg: %.4f\n', results.auprc_neg, results.auroc_neg);
% fprintf('Pearson r (pos edges): %.4f   (neg edges): %.4f   (overall): %.4f\n', ...
%     results.r_pos, results.r_neg, results.r_overall);
% fprintf('MAE (overall): %.4f  MAE (pos true edges): %.4f  MAE (neg true edges): %.4f\n', ...
%     results.mae_overall, results.mae_pos, results.mae_neg);
% fprintf('Sign accuracy (on true edges): pos %.3f  neg %.3f  overall %.3f\n', ...
%     results.sign_acc_pos, results.sign_acc_neg, results.sign_acc_overall);
% fprintf('Pred density: %.4f  True density: %.4f\n', results.pred_density, results.true_density);
% fprintf('SHD (adjacency mismatch): %d   Signed-SHD (sign mismatch counted): %d\n', ...
%     results.SHD, results.SHD_signed);

%% ----------------- local function -----------------
function out = compare_connectivity_matrices(A, B, nboot)
% returns struct with many metrics comparing A (ground-truth) and B (predicted)
% A, B: NxN numeric matrices (directed, can be signed)
% nboot: number of bootstrap resamples over off-diagonal entries for CIs (0 to skip)

N = size(A,1);
offdiag = ~eye(N);

% vectors of off-diagonal entries
vecA = A(offdiag);
vecB = B(offdiag);

% true labels
labels_pos = double(vecA > 0);   % 1 for excitatory edges
labels_neg = double(vecA < 0);   % 1 for inhibitory edges

% scores for ranking (higher => more likely excitatory). For inhibitory tests we'll flip sign.
scores = vecB;

% AUPRC / AUROC (positive)
if numel(unique(labels_pos)) > 1
    [~,~,~,auroc_pos] = perfcurve(labels_pos, scores, 1);
    [~,~,~,auprc_pos] = perfcurve(labels_pos, scores, 1, 'xCrit','reca','yCrit','prec');
else
    auroc_pos = NaN; auprc_pos = NaN;
end

% AUPRC / AUROC (negative) — treat negative edges as positive class after flipping scores
if numel(unique(labels_neg)) > 1
    [~,~,~,auroc_neg] = perfcurve(labels_neg, -scores, 1);
    [~,~,~,auprc_neg] = perfcurve(labels_neg, -scores, 1, 'xCrit','reca','yCrit','prec');
else
    auroc_neg = NaN; auprc_neg = NaN;
end

% Pearson r on true edges (pos / neg)
mask_pos = (A > 0) & offdiag;
mask_neg = (A < 0) & offdiag;

if any(mask_pos(:))
    r_pos = corr(B(mask_pos), A(mask_pos));
    mae_pos = mean(abs(B(mask_pos) - A(mask_pos)));
    rmse_pos = sqrt(mean((B(mask_pos) - A(mask_pos)).^2));
    sign_acc_pos = mean(sign(B(mask_pos)) == sign(A(mask_pos)));
else
    r_pos = NaN; mae_pos = NaN; rmse_pos = NaN; sign_acc_pos = NaN;
end

if any(mask_neg(:))
    r_neg = corr(B(mask_neg), A(mask_neg));
    mae_neg = mean(abs(B(mask_neg) - A(mask_neg)));
    rmse_neg = sqrt(mean((B(mask_neg) - A(mask_neg)).^2));
    sign_acc_neg = mean(sign(B(mask_neg)) == sign(A(mask_neg)));
else
    r_neg = NaN; mae_neg = NaN; rmse_neg = NaN; sign_acc_neg = NaN;
end

% overall correlations and errors (off-diagonal)
r_overall = corr(vecA, vecB, 'rows','complete');
mae_overall = mean(abs(vecB - vecA));
rmse_overall = sqrt(mean((vecB - vecA).^2));
sign_acc_overall = mean(sign(vecB) == sign(vecA)); % counts zeros as equal-sign if both zero

% densities (off-diag)
true_density = mean(vecA ~= 0);
pred_density = mean(vecB ~= 0);

% predicted top-K matching: pick top-Kpos excitatory and top-Kneg inhibitory
Kpos = sum(labels_pos == 1);
Kneg = sum(labels_neg == 1);

pred_adj = false(N); % predicted adjacency (binary)
% create indices sorted by score for off-diagonal linear indexing
[s_sorted_desc, idx_desc] = sort(scores, 'descend');
[s_sorted_asc, idx_asc] = sort(scores, 'ascend');

if Kpos > 0
    sel_pos = idx_desc(1:min(Kpos, numel(idx_desc)));
else
    sel_pos = [];
end
if Kneg > 0
    sel_neg = idx_asc(1:min(Kneg, numel(idx_asc)));
else
    sel_neg = [];
end

% map linear off-diagonal indices back to full NxN indices
[rows_off, cols_off] = find(offdiag);
% sel_pos and sel_neg are indices into vec arrays -> map:
pred_adj(sub2ind([N,N], rows_off(sel_pos), cols_off(sel_pos))) = true;
pred_adj(sub2ind([N,N], rows_off(sel_neg), cols_off(sel_neg))) = true;

% true adjacency (binary)
true_adj = (A ~= 0);

% SHD: number of differing edges (on off-diagonal)
SHD = sum(xor(true_adj(offdiag), pred_adj(offdiag)));

% signed SHD: count where adjacency is same but sign differs, plus adjacency mismatches
% We'll compute as adjacency mismatch (in/out) + sign mismatches on common edges
adj_mismatch = sum(xor(true_adj(offdiag), pred_adj(offdiag)));
common_mask = true_adj & pred_adj; % edges predicted and truly present
signed_mismatch = 0;
if any(common_mask(:))
    signed_mismatch = sum(sign(B(common_mask)) ~= sign(A(common_mask)));
end
SHD_signed = adj_mismatch + signed_mismatch;

% precision/recall/F1 at matching-density threshold (combined)
pred_binary_vec = pred_adj(offdiag);
true_binary_vec = true_adj(offdiag);

tp = sum(pred_binary_vec & true_binary_vec);
fp = sum(pred_binary_vec & ~true_binary_vec);
fn = sum(~pred_binary_vec & true_binary_vec);

precision = tp / (tp + fp + eps);
recall = tp / (tp + fn + eps);
F1 = 2 * precision * recall / (precision + recall + eps);

% weight correlations on sign-separated true edges (already r_pos, r_neg)
% also compute correlation on predicted edges
if any(pred_adj(offdiag))
    r_pred_on_predicted_edges = corr(vecA(pred_adj(offdiag)), vecB(pred_adj(offdiag)), 'rows','complete');
else
    r_pred_on_predicted_edges = NaN;
end

% Bootstrapped CI for auprc_pos and auprc_neg if requested
if nboot > 0
    rng(0); % reproducible
    n_off = numel(vecA);
    auprc_pos_bs = nan(nboot,1);
    auprc_neg_bs = nan(nboot,1);
    for b = 1:nboot
        idx = randsample(n_off, n_off, true);
        lb = labels_pos(idx);
        sb = scores(idx);
        if numel(unique(lb)) > 1
            [~,~,~,auprc_pos_bs(b)] = perfcurve(lb, sb, 1, 'xCrit','reca','yCrit','prec');
        else
            auprc_pos_bs(b) = NaN;
        end
        ln = labels_neg(idx);
        if numel(unique(ln)) > 1
            [~,~,~,auprc_neg_bs(b)] = perfcurve(ln, -sb, 1, 'xCrit','reca','yCrit','prec');
        else
            auprc_neg_bs(b) = NaN;
        end
    end
    auprc_pos_CI = prctile(auprc_pos_bs, [2.5 97.5]);
    auprc_neg_CI = prctile(auprc_neg_bs, [2.5 97.5]);
else
    auprc_pos_CI = [NaN NaN];
    auprc_neg_CI = [NaN NaN];
end

% pack outputs
out.auprc_pos = auprc_pos;
out.auprc_neg = auprc_neg;
out.auroc_pos = auroc_pos;
out.auroc_neg = auroc_neg;
out.auprc_pos_CI = auprc_pos_CI;
out.auprc_neg_CI = auprc_neg_CI;

out.r_pos = r_pos;
out.r_neg = r_neg;
out.r_overall = r_overall;
out.r_pred_on_predicted_edges = r_pred_on_predicted_edges;

out.mae_pos = mae_pos;
out.mae_neg = mae_neg;
out.mae_overall = mae_overall;
out.rmse_pos = rmse_pos;
out.rmse_neg = rmse_neg;
out.rmse_overall = rmse_overall;

out.sign_acc_pos = sign_acc_pos;
out.sign_acc_neg = sign_acc_neg;
out.sign_acc_overall = sign_acc_overall;

out.true_density = true_density;
out.pred_density = pred_density;

out.SHD = SHD;
out.SHD_signed = SHD_signed;

out.precision_at_match = precision;
out.recall_at_match = recall;
out.F1_at_match = F1;

% create a small table for easy display
varnames = {'auprc_pos','auroc_pos','auprc_neg','auroc_neg', ...
            'r_pos','r_neg','r_overall','mae_overall','rmse_overall','sign_acc_overall', ...
            'true_density','pred_density','SHD','SHD_signed','precision_match','recall_match','F1_match'};
vartable = {out.auprc_pos, out.auroc_pos, out.auprc_neg, out.auroc_neg, ...
            out.r_pos, out.r_neg, out.r_overall, out.mae_overall, out.rmse_overall, out.sign_acc_overall, ...
            out.true_density, out.pred_density, out.SHD, out.SHD_signed, out.precision_at_match, out.recall_at_match, out.F1_at_match};

out.summary_table = cell2table(vartable, 'VariableNames', varnames);

end