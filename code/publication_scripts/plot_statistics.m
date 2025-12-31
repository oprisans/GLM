clear all; close all; clc;
warning('off');

% ---------------- I/O and paths ----------------
fullpath = mfilename('fullpath');
if isempty(fullpath)
    fullpath = matlab.desktop.editor.getActiveFilename;
end
[scriptFolder, scriptName, ~] = fileparts(fullpath);
parentFolder = fileparts(scriptFolder);

dataPath   = fullfile(parentFolder, 'results','zScore_seed_1_50');
resultsPath= fullfile(parentFolder, 'results','zScore_seed_1_50');
if ~exist(resultsPath, 'dir'), mkdir(resultsPath); end

%fileName=fullfile(dataPath,'a1_univariate_1y_1x_lassoglm_groundTruth_zScore_comp_groundTruth.xlsx');
fileName=fullfile(dataPath,'aa1_lassoglm_1y_ALLx_zScore_comp_groundTruth_seed_1_50.xlsx');

% 
opts = detectImportOptions(fileName,'Sheet','b1_CV1');
opts.VariableNamingRule = 'preserve';      % keep column names exactly as in Excel
% (optionally) adjust types:
% opts = setvartype(opts,{'N','T','sparsity'},'double');

T = readtable(fileName,opts);
% fixed values
% INPUT: T (table loaded with readtable)
plot_paper_figures(T, resultsPath)


%%%%%%%%%%%%%%
function plot_paper_figures(T, saveDir)
% plot_paper_figures  Produce paper-ready figures from results table
% Usage:
%   plot_paper_figures(T)            % T is a MATLAB table
%   plot_paper_figures()             % will try to read 'results.xlsx'
%   plot_paper_figures(T,'figs')     % save figures to ./figs
%
% Output: four figures (AUPRC vs sparsity; F1/precision/recall vs sparsity;
%        AUPRC & F1 vs snr_db; r_overall & mae_overall vs T)
%
% Notes:
% - Requires columns named as in your listing (auprc_pos, F1_at_match, etc).
% - Uses bootstrap (1000) to compute 95% CI when per-group sample size >=2.
% - Change Nlist at top if you want different Ns.

% -----------------------
% Parameters (change here if desired)
Nlist = [10 20 30];
% fixed parameter values
fixed.T = 2400;
fixed.dt = 0.05;
fixed.snr_db = 0;
fixed.sigma_r = 0.2;
fixed.tau = 1.2;
fixed.w = 0.25;
fixed.rho_target = 0.9;
fixed.smooth_sigma = 0.2;
fixed.seed = 42;
tol = 1e-6;
nboot = 1000;
% -----------------------

if nargin < 1 || isempty(T)
    % try to read from default filename
    fname = 'results.xlsx';
    if ~isfile(fname)
        error('No table T provided and %s not found in current folder.', fname);
    end
    T = readtable(fname);
end

if nargin < 2
    saveDir = '';
end
if ~isempty(saveDir) && ~isfolder(saveDir)
    mkdir(saveDir);
end

% check required columns
requiredCols = {'auprc_pos','auprc_pos_CI_low','auprc_pos_CI_high', ...
    'F1_at_match','precision_at_match','recall_at_match', ...
    'r_overall','mae_overall','true_density','sparsity','snr_db','N','T'};
missing = setdiff(requiredCols, T.Properties.VariableNames);
% Not fatal if CI columns missing; allow bootstrap fallback
missingRequired = setdiff({'auprc_pos','F1_at_match','r_overall','mae_overall','sparsity','snr_db','N','T'}, T.Properties.VariableNames);
if ~isempty(missingRequired)
    error('Table is missing required columns: %s', strjoin(missingRequired,', '));
end

% Helper function: compute mean and CI per unique x
    function [xvals, mu, ci_lo, ci_hi, n] = group_stats(tbl, xvar, statName)
        % handles numeric or categorical xvar
        xcol = tbl.(xvar);
        xvals = unique(xcol);
        % convert to numeric for sorting if possible
        if isnumeric(xvals) || islogical(xvals)
            [xvals, ixsort] = sort(double(xvals));
        else
            xvals = unique(string(xcol),'stable');
            ixsort = 1:numel(xvals);
        end
        xvals = xvals(ixsort);
        m = numel(xvals);
        mu = nan(m,1); ci_lo = nan(m,1); ci_hi = nan(m,1); n = zeros(m,1);
        for i = 1:m
            xv = xvals(i);
            if isnumeric(tbl.(xvar))
                sel = abs(tbl.(xvar) - double(xv)) <= tol;
            else
                sel = string(tbl.(xvar)) == string(xv);
            end
            vals = tbl{sel, statName};
            vals = vals(~isnan(vals));
            n(i) = numel(vals);
            if n(i)==0
                mu(i)=NaN; ci_lo(i)=NaN; ci_hi(i)=NaN;
                continue;
            end
            mu(i) = mean(vals);
            if n(i) >= 2
                % bootstrap mean
                b = bootstrp(nboot, @mean, vals);
                ci = prctile(b, [2.5 97.5]);
                ci_lo(i) = ci(1);
                ci_hi(i) = ci(2);
            else
                % fallback to SEM (n == 1)
                sem = std(vals)/sqrt(max(1,n(i)));
                ci_lo(i) = mu(i) - 1.96*sem;
                ci_hi(i) = mu(i) + 1.96*sem;
            end
        end
    end

% Helper: shaded CI plotting (asymmetric)
    function plot_shaded(x, mu, lo, hi, lineSpec)
        % remove NaNs for patch
        valid = ~isnan(mu) & ~isnan(lo) & ~isnan(hi);
        if all(~valid), return; end
        x = double(x(valid)); mu = mu(valid); lo = lo(valid); hi = hi(valid);
        xx = [x; flipud(x)];
        yy = [lo; flipud(hi)];
        p = patch(xx, yy, 1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % set color to current axes color order
        color = get(gca,'ColorOrder'); color = color(mod(get(gca,'NextPlot')-2,size(color,1))+1,:);
        set(p,'FaceColor',color(1,:));
        hold on;
        plot(x, mu, lineSpec, 'LineWidth', 1.5, 'MarkerSize', 6);
    end

% -----------------------
% FIGURE 1: AUPRC_pos vs sparsity (overlay N)
% -----------------------
figure('Name','Fig1_AUPRC_vs_sparsity','NumberTitle','off');
hold on;
legendEntries = {};
for k = 1:numel(Nlist)
    Nval = Nlist(k);
    % mask fixed params except sparsity
    mask = abs(T.N - Nval) <= tol ...
        & abs(T.T - fixed.T) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.snr_db - fixed.snr_db) <= tol ...
        & abs(T.sigma_r - fixed.sigma_r) <= tol ...
        & abs(T.tau - fixed.tau) <= tol ...
        & abs(T.w - fixed.w) <= tol ...
        & abs(T.rho_target - fixed.rho_target) <= tol ...
        & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
        & abs(T.seed - fixed.seed) <= tol;
    rows = T(mask,:);
    if isempty(rows)
        fprintf('Fig1: no rows for N=%d. skipping\n', Nval);
        continue;
    end
    [xvals, mu, ci_lo, ci_hi, n] = group_stats(rows, 'sparsity', 'auprc_pos');
    % If explicit CI columns exist, prefer their aggregated mean of endpoints
    if ismember('auprc_pos_CI_low', T.Properties.VariableNames) && ismember('auprc_pos_CI_high', T.Properties.VariableNames)
        % compute per-group mean CI endpoints
        m = numel(xvals);
        for i=1:m
            xv = xvals(i);
            if isnumeric(rows.sparsity)
                sel = abs(rows.sparsity - double(xv)) <= tol;
            else
                sel = string(rows.sparsity) == string(xv);
            end
            if any(sel)
                lo_col = rows.auprc_pos_CI_low(sel); hi_col = rows.auprc_pos_CI_high(sel);
                lo_col = lo_col(~isnan(lo_col)); hi_col = hi_col(~isnan(hi_col));
                if ~isempty(lo_col) && ~isempty(hi_col)
                    mu(i) = mean(rows.auprc_pos(sel)); % ensure mean
                    ci_lo(i) = mean(lo_col);
                    ci_hi(i) = mean(hi_col);
                end
            end
        end
    end
    % plot shaded CI + mean
    plot_shaded(xvals, mu, ci_lo, ci_hi, '-o');
    legendEntries{end+1} = sprintf('N=%d (n_{avg}=%.1f)', Nval, mean(n)); %#ok<SAGROW>
end
% plot prevalence baseline if present (true_density)
% try to compute prevalence per sparsity from T (average across Ns)
try
    % pick reference Nlist(1) with same mask as earlier
    refN = Nlist(1);
    refMask = abs(T.N - refN) <= tol ...
        & abs(T.T - fixed.T) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.snr_db - fixed.snr_db) <= tol;
    refRows = T(refMask,:);
    if ~isempty(refRows)
        [svals, ~, ~] = unique(refRows.sparsity);
        % for each sparsity compute average true_density
        baseline = nan(size(svals));
        for i=1:numel(svals)
            if isnumeric(refRows.sparsity)
                sel = abs(refRows.sparsity - double(svals(i))) <= tol;
            else
                sel = string(refRows.sparsity) == string(svals(i));
            end
            baseline(i) = mean(refRows.true_density(sel));
        end
        plot(double(svals), baseline, '--k', 'LineWidth', 1);
        legendEntries{end+1} = 'prevalence (baseline)';
    end
end
xlabel('sparsity'); ylabel('AUPRC_{pos}'); title('AUPRC_{pos} vs sparsity (95% CI)');
grid on; legend(legendEntries,'Location','best'); hold off;
if ~isempty(saveDir), saveas(gcf, fullfile(saveDir,'Fig1_AUPRC_vs_sparsity.png')); end

% -----------------------
% FIGURE 2: F1_at_match vs sparsity + precision/recall panel
% -----------------------
figure('Name','Fig2_F1_precision_recall','NumberTitle','off');
% left: F1 vs sparsity
subplot(1,2,1); hold on;
legendEntries = {};
for k = 1:numel(Nlist)
    Nval = Nlist(k);
    mask = abs(T.N - Nval) <= tol ...
        & abs(T.T - fixed.T) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.snr_db - fixed.snr_db) <= tol ...
        & abs(T.sigma_r - fixed.sigma_r) <= tol ...
        & abs(T.tau - fixed.tau) <= tol ...
        & abs(T.w - fixed.w) <= tol ...
        & abs(T.rho_target - fixed.rho_target) <= tol ...
        & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
        & abs(T.seed - fixed.seed) <= tol;
    rows = T(mask,:);
    if isempty(rows), continue; end
    [xvals, mu, ci_lo, ci_hi, n] = group_stats(rows, 'sparsity', 'F1_at_match');
    plot_shaded(xvals, mu, ci_lo, ci_hi, '-s');
    legendEntries{end+1} = sprintf('N=%d', Nval); %#ok<SAGROW>
end
xlabel('sparsity'); ylabel('F1_{at\_match}'); title('F1_{at\_match} vs sparsity'); grid on; legend(legendEntries,'Location','best'); hold off;

% right: precision & recall (pick first N in Nlist as representative)
subplot(1,2,2); hold on;
repN = Nlist(1);
mask = abs(T.N - repN) <= tol ...
    & abs(T.T - fixed.T) <= tol ...
    & abs(T.dt - fixed.dt) <= tol ...
    & abs(T.sigma_r - fixed.sigma_r) <= tol ...
    & abs(T.tau - fixed.tau) <= tol ...
    & abs(T.w - fixed.w) <= tol ...
    & abs(T.rho_target - fixed.rho_target) <= tol ...
    & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
    & abs(T.seed - fixed.seed) <= tol;
rows = T(mask,:);
if ~isempty(rows)
    [xvals_p, mu_p, lo_p, hi_p, n_p] = group_stats(rows, 'sparsity', 'precision_at_match');
    [xvals_r, mu_r, lo_r, hi_r, n_r] = group_stats(rows, 'sparsity', 'recall_at_match');
    plot_shaded(xvals_p, mu_p, lo_p, hi_p, '-o');
    plot_shaded(xvals_r, mu_r, lo_r, hi_r, '-^');
    legend({sprintf('precision (N=%d)',repN), sprintf('recall (N=%d)',repN)}, 'Location','best');
end
xlabel('sparsity'); ylabel('score'); title(sprintf('precision/recall (N=%d)', repN)); grid on; hold off;
if ~isempty(saveDir), saveas(gcf, fullfile(saveDir,'Fig2_F1_precision_recall.png')); end

% -----------------------
% FIGURE 3: noise sensitivity (AUPRC_pos and F1_at_match vs snr_db)
% fixed sparsity pick typical value (use median of available sparsities)
% -----------------------
% pick a sparsity value present in the data (median)
allSpars = unique(T.sparsity);
if isempty(allSpars)
    sparsity_fixed = 0.05;
else
    sparsity_fixed = median(allSpars);
end

figure('Name','Fig3_noise_sensitivity','NumberTitle','off');
subplot(2,1,1); hold on;
legendEntries = {};
for k = 1:numel(Nlist)
    Nval = Nlist(k);
    mask = abs(T.N - Nval) <= tol ...
        & abs(T.T - fixed.T) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.sigma_r - fixed.sigma_r) <= tol ...
        & abs(T.tau - fixed.tau) <= tol ...
        & abs(T.w - fixed.w) <= tol ...
        & abs(T.rho_target - fixed.rho_target) <= tol ...
        & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
        & abs(T.seed - fixed.seed) <= tol ...
        & abs(T.sparsity - sparsity_fixed) <= tol;
    rows = T(mask,:);
    if isempty(rows), continue; end
    [xvals, mu, lo, hi, n] = group_stats(rows, 'snr_db', 'auprc_pos');
    plot_shaded(xvals, mu, lo, hi, '-o');
    legendEntries{end+1} = sprintf('N=%d',Nval); %#ok<SAGROW>
end
xlabel('snr\_db'); ylabel('AUPRC_{pos}'); title(sprintf('AUPRC_{pos} vs SNR (sparsity=%.3g)', sparsity_fixed)); grid on; legend(legendEntries,'Location','best'); hold off;

subplot(2,1,2); hold on;
legendEntries = {};
for k = 1:numel(Nlist)
    Nval = Nlist(k);
    mask = abs(T.N - Nval) <= tol ...
        & abs(T.T - fixed.T) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.sigma_r - fixed.sigma_r) <= tol ...
        & abs(T.tau - fixed.tau) <= tol ...
        & abs(T.w - fixed.w) <= tol ...
        & abs(T.rho_target - fixed.rho_target) <= tol ...
        & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
        & abs(T.seed - fixed.seed) <= tol ...
        & abs(T.sparsity - sparsity_fixed) <= tol;
    rows = T(mask,:);
    if isempty(rows), continue; end
    [xvals, mu, lo, hi, n] = group_stats(rows, 'snr_db', 'F1_at_match');
    plot_shaded(xvals, mu, lo, hi, '-s');
    legendEntries{end+1} = sprintf('N=%d',Nval); %#ok<SAGROW>
end
xlabel('snr\_db'); ylabel('F1_{at\_match}'); title(sprintf('F1_{at\_match} vs SNR (sparsity=%.3g)', sparsity_fixed)); grid on; legend(legendEntries,'Location','best'); hold off;
if ~isempty(saveDir), saveas(gcf, fullfile(saveDir,'Fig3_noise_sensitivity.png')); end

% -----------------------
% FIGURE 4: r_overall and mae_overall vs T (log x-axis)
% -----------------------
figure('Name','Fig4_data_requirements','NumberTitle','off');
subplot(2,1,1); hold on;
legendEntries = {};
for k = 1:numel(Nlist)
    Nval = Nlist(k);
    mask = abs(T.N - Nval) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.snr_db - fixed.snr_db) <= tol ...
        & abs(T.sigma_r - fixed.sigma_r) <= tol ...
        & abs(T.tau - fixed.tau) <= tol ...
        & abs(T.w - fixed.w) <= tol ...
        & abs(T.rho_target - fixed.rho_target) <= tol ...
        & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
        & abs(T.seed - fixed.seed) <= tol;
    rows = T(mask,:);
    if isempty(rows), continue; end
    [xvals, mu, lo, hi, n] = group_stats(rows, 'T', 'r_overall');
    % sort by xvals
    [xvals, sidx] = sort(double(xvals)); mu = mu(sidx); lo = lo(sidx); hi = hi(sidx);
    plot_shaded(xvals, mu, lo, hi, '-o');
    legendEntries{end+1} = sprintf('N=%d',Nval); %#ok<SAGROW>
end
set(gca,'XScale','log');
xlabel('T (log scale)'); ylabel('r\_overall'); title('r\_overall vs T (data length)'); grid on; legend(legendEntries,'Location','best'); hold off;

subplot(2,1,2); hold on;
legendEntries = {};
for k = 1:numel(Nlist)
    Nval = Nlist(k);
    mask = abs(T.N - Nval) <= tol ...
        & abs(T.dt - fixed.dt) <= tol ...
        & abs(T.snr_db - fixed.snr_db) <= tol ...
        & abs(T.sigma_r - fixed.sigma_r) <= tol ...
        & abs(T.tau - fixed.tau) <= tol ...
        & abs(T.w - fixed.w) <= tol ...
        & abs(T.rho_target - fixed.rho_target) <= tol ...
        & abs(T.smooth_sigma - fixed.smooth_sigma) <= tol ...
        & abs(T.seed - fixed.seed) <= tol;
    rows = T(mask,:);
    if isempty(rows), continue; end
    [xvals, mu, lo, hi, n] = group_stats(rows, 'T', 'mae_overall');
    [xvals, sidx] = sort(double(xvals)); mu = mu(sidx); lo = lo(sidx); hi = hi(sidx);
    plot_shaded(xvals, mu, lo, hi, '-s');
    legendEntries{end+1} = sprintf('N=%d',Nval); %#ok<SAGROW>
end
set(gca,'XScale','log');
xlabel('T (log scale)'); ylabel('MAE_{overall}'); title('MAE_{overall} vs T'); grid on; legend(legendEntries,'Location','best'); hold off;
if ~isempty(saveDir), saveas(gcf, fullfile(saveDir,'Fig4_data_requirements.png')); end

% End
fprintf('Figures generated. Save dir: %s\n', saveDir);

end