% compare_groundTruth_vs_aa1_lassoglm_1x_ALLy_zScore.m
% Compare predicted coefficient outputs from aa1_lassoglm against ground-truth.
% Parallelized over metric files produced by aa1_lassoglm_1y_ALLx_zScore.
%
% Optional in base workspace before calling:
%   seedList = [42 101];    % process only these seeds
%   seedval  = 42;          % single seed
%   nboot = 200;            % bootstrap count (default 200)
%
% Writes per-file compare_*.xlsx in results/zScore and aggregated compare_*.xlsx/.mat

clearvars -except seedList seedval nboot
if ~exist('nboot','var'), nboot = 200; end
tol = 1e-6;

% --- resolve folders relative to script location ---
parentFolder = fileparts(fileparts(mfilename('fullpath')));
if isempty(parentFolder), parentFolder = pwd; end
dataPath = fullfile(parentFolder, 'data');
outResultsDir = fullfile(parentFolder, 'results', 'zScore');
if ~exist(outResultsDir,'dir'), mkdir(outResultsDir); end

% metrics produced by aa1 (expected location)
metricsDir = fullfile(parentFolder, 'data', 'aa1_lassoglm_1y_ALLx_zScore');
if ~isfolder(metricsDir)
    % fallback to data/<scriptname> if different naming used
    metricsDir = fullfile(parentFolder, 'data', 'aa1_lassoglm_1y_1x_zScore');
end
if ~isfolder(metricsDir)
    warning('Metrics directory not found: %s. Aborting.', metricsDir);
    return;
end

% index ground-truth files (xlsx or mat) under data/groundTruth
gt_files = {};
gtDir = fullfile(parentFolder, 'data', 'groundTruth');
if isfolder(gtDir)
    fx = dir(fullfile(gtDir,'**','*.xlsx'));
    fm = dir(fullfile(gtDir,'**','*.mat'));
    for k=1:numel(fx), gt_files{end+1} = fullfile(fx(k).folder, fx(k).name); end
    for k=1:numel(fm), gt_files{end+1} = fullfile(fm(k).folder, fm(k).name); end
end
fprintf('Indexed %d ground-truth files.\n', numel(gt_files));

% discover metric files (xlsx and mat)
mx = dir(fullfile(metricsDir,'aa1_lassoglm_1y_ALLx_zScore_*.xlsx'));
mm = dir(fullfile(metricsDir,'aa1_lassoglm_1y_ALLx_zScore_*.mat'));
metric_files = [mx; mm];
if isempty(metric_files)
    % try more permissive pattern
    mx = dir(fullfile(metricsDir,'*.xlsx'));
    mm = dir(fullfile(metricsDir,'*.mat'));
    metric_files = [mx; mm];
end
metric_paths = fullfile({metric_files.folder}, {metric_files.name})';

% optional seed filter from base
seedFilter = [];
if exist('seedList','var') && ~isempty(seedList)
    seedFilter = double(seedList(:))';
elseif exist('seedval','var') && ~isempty(seedval)
    seedFilter = double(seedval);
end
if ~isempty(seedFilter)
    keep = false(size(metric_paths));
    for i=1:numel(metric_paths)
        fn = metric_files(i).name;
        tok = regexp(fn,'seed[_-]?([0-9]+)','tokens','once');
        if ~isempty(tok) && any(str2double(tok{1})==seedFilter)
            keep(i) = true;
        end
    end
    metric_paths = metric_paths(keep);
end

if isempty(metric_paths)
    warning('No metric files found after filtering. Exiting.');
    return;
end
fprintf('Found %d metric files to compare.\n', numel(metric_paths));

% candidate sheet names / variable names for coefficient matrices
candidateSheets = {'b1_CV1','b1_CV2','b1_AIC','b1_BIC','b1_Custom','b1_custom','b1CV1','b1CV2','b1_Aic','b1_Bic'};

% Process files in parallel. Each iteration writes a per-file compare_<base>_comp_groundTruth.xlsx
pool = gcp('nocreate');
if ~isempty(pool)
    pctRunOnAll addpath(fileparts(mfilename('fullpath')));
end

parfor fi = 1:numel(metric_paths)
    metricPath = metric_paths{fi};
    [~, base, ext] = fileparts(metricPath);
    try
        fprintf('Worker %d processing %s\n', getWorkerID(), base);
    catch
        fprintf('Processing %s\n', base);
    end

    % extract keys from metric filename (N_, T_, dt_, sparsity_, snr_, seed_)
    fname = lower([base, ext]);
    tokN = regexp(fname,'n[_-]?([0-9]+)','tokens','once');
    tokT = regexp(fname,'t[_-]?([0-9]+)','tokens','once');
    tokdt = regexp(fname,'dt[_-]?([0-9\.]+)','tokens','once');
    tokspars = regexp(fname,'sparsity[_-]?([0-9\.]+)','tokens','once');
    toksnr = regexp(fname,'snr[_-]?([0-9]+)','tokens','once');
    tokseed = regexp(fname,'seed[_-]?([0-9]+)','tokens','once');

    % build keys for ground-truth matching, use flexible formatting
    keys = {};
    if ~isempty(tokN), keys{end+1} = sprintf('n_%s', tokN{1}); end
    if ~isempty(tokT), keys{end+1} = sprintf('t_%s', tokT{1}); end
    if ~isempty(tokdt), keys{end+1} = sprintf('dt_%s', strrep(tokdt{1},'.',',') ); end
    if ~isempty(tokspars), keys{end+1} = sprintf('sparsity_%s', strrep(tokspars{1},'.',',') ); end
    if ~isempty(toksnr), keys{end+1} = sprintf('snr_%s', toksnr{1}); end
    if ~isempty(tokseed), keys{end+1} = sprintf('seed_%s', tokseed{1}); end

    % find ground-truth file matching keys
    gt_match = find_file_by_keys(gt_files, keys);
    if isempty(gt_match)
        warning('No ground-truth match for metric %s. Skipping.', base);
        continue;
    end

    % read ground-truth A matrix
    [A_true, srcA] = robust_read_groundtruth(gt_match);
    if isempty(A_true)
        warning('Failed to read ground-truth for %s. Skipping.', base);
        continue;
    end

    % read predicted matrices (try candidate sheets/vars)
    Bmats = struct();
    for s = 1:numel(candidateSheets)
        try
            tmp = robust_read_pred(metricPath, {candidateSheets{s}});
            if ~isempty(tmp)
                % normalize name to a tag key
                tag = candidateSheets{s};
                tag = regexprep(tag,'[^a-zA-Z0-9]','_');
                Bmats.(tag) = tmp;
            end
        catch
        end
    end

    % if nothing read, skip
    fld = fieldnames(Bmats);
    if isempty(fld)
        warning('No predicted coefficient sheet found in %s. Skipping.', base);
        continue;
    end

    % for each matrix compute metrics
    rows = [];
    methods = {};
    for kf = 1:numel(fld)
        method = fld{kf};
        Braw = Bmats.(method);
        try
            Bmat = ensure_square_NxN(Braw, size(A_true,1));
        catch
            warning('Cannot normalize predicted matrix for %s method %s. Skipping this method.', base, method);
            continue;
        end
        out = compare_connectivity_matrices(A_true, Bmat, nboot, tol);

        % attempt to parse params from keys (fallback to NaN)
        Nval = try_parse_key_number(keys, 'n_'); 
        Tval = try_parse_key_number(keys, 't_');
        dtval = try_parse_key_number(keys, 'dt_');
        sparsityval = try_parse_key_number(keys, 'sparsity_');
        snr_dbval = try_parse_key_number(keys, 'snr_');
        seedval_local = try_parse_key_number(keys, 'seed_');

        % pack row
        row = pack_row(Nval, Tval, dtval, sparsityval, snr_dbval, seedval_local, ...
                       nan, nan, nan, nan, nan, out); % w,sigma,tau,rho,smooth are NaN here
        rows = [rows; row]; %#ok<AGROW>
        methods{end+1,1} = method; %#ok<AGROW>
    end

    if isempty(rows)
        continue;
    end

    % make table and add method column
    variableNames = {'N','T','dt','sparsity','snr_db','seed','w','sigma_r','tau','rho_target','smooth_sigma', ...
                     'auprc_pos','auprc_neg','auroc_pos','auroc_neg','auprc_pos_CI_low','auprc_pos_CI_high', ...
                     'auprc_neg_CI_low','auprc_neg_CI_high','r_pos','r_neg','r_overall','r_pred_on_predicted_edges', ...
                     'mae_pos','mae_neg','mae_overall','rmse_pos','rmse_neg','rmse_overall','sign_acc_pos','sign_acc_neg', ...
                     'sign_acc_overall','true_density','pred_density','SHD','SHD_signed','precision_at_match','recall_at_match','F1_at_match'};
    T = array2table(rows, 'VariableNames', variableNames);
    T.method = methods;

    % write per-file compare output
    outBase = sprintf('compare_%s_comp_groundTruth.xlsx', base);
    outPath = fullfile(outResultsDir, outBase);
    try
        writetable(T, outPath, 'Sheet', 'metrics');
    catch ME
        warning('Failed to write per-file compare for %s: %s', base, ME.message);
    end
end % parfor

% ---------------- collect per-file compare tables and aggregate per-method -----------
cmp_files = dir(fullfile(outResultsDir,'compare_*_comp_groundTruth.xlsx'));
if isempty(cmp_files)
    fprintf('No per-file compare outputs found in %s\n', outResultsDir);
    return;
end

allT = table();
for k = 1:numel(cmp_files)
    try
        tt = readtable(fullfile(cmp_files(k).folder, cmp_files(k).name), 'Sheet', 'metrics');
        if ~isempty(tt)
            allT = [allT; tt]; %#ok<AGROW>
        end
    catch
        warning('Failed reading %s', cmp_files(k).name);
    end
end

if isempty(allT)
    warning('No readable compare tables. Exiting.');
    return;
end

% normalize names to match old pipeline expectations
allT.Properties.VariableNames = matlab.lang.makeValidName(allT.Properties.VariableNames);

% split by method and write aggregated files similar to original output
methods = unique(allT.method);
for m = 1:numel(methods)
    meth = methods{m};
    sel = strcmp(allT.method, meth);
    Tm = allT(sel, :);
    if isempty(Tm), continue; end
    % remove method column for compatibility and reorder columns to original variableNames
    if ismember('method', Tm.Properties.VariableNames)
        Tm.method = [];
    end
    % ensure columns order matches variableNames if present
    varsToKeep = intersect(variableNames, Tm.Properties.VariableNames, 'stable');
    Tm = Tm(:, varsToKeep);

    tag = meth; % e.g. b1_CV1 etc.
    % sanitize tag for filenames
    tagSafe = regexprep(tag,'[^A-Za-z0-9_]','_');
    fname_xlsx = fullfile(outResultsDir, sprintf('compare_lassoglm_b1_%s_comp_groundTruth.xlsx', tagSafe));
    saveMat = fullfile(outResultsDir, sprintf('compare_lassoglm_b1_%s_comp_groundTruth.mat', tagSafe));
    try
        writetable(Tm, fname_xlsx, 'Sheet', tagSafe);
        T = Tm; %#ok<NASGU>
        save(saveMat, 'T');
        fprintf('Wrote aggregated %s rows: %d -> %s\n', tagSafe, height(Tm), fname_xlsx);
    catch ME
        warning('Failed writing aggregated %s: %s', tagSafe, ME.message);
    end
end

fprintf('Done. Results (if any) are in %s\n', outResultsDir);

%% ---------------- helpers (reused from your prior script) ----------------

function f = find_file_by_keys(fileList, keys)
    f = '';
    if isempty(fileList), return; end
    if isempty(keys), return; end
    kl = lower(keys);
    for i=1:numel(fileList)
        nm = lower(fileList{i});
        match = true;
        for j=1:numel(kl)
            if ~contains(nm, kl{j})
                match = false; break;
            end
        end
        if match, f = fileList{i}; return; end
    end
end

function [M, src] = robust_read_groundtruth(pathOrFile)
    M = []; src = '';
    if isempty(pathOrFile), return; end
    [~,~,ext] = fileparts(pathOrFile);
    try
        if strcmpi(ext,'.xlsx') || strcmpi(ext,'.xls')
            sheetsToTry = {'A','A_true','Adj','adj','Sheet1','a'};
            for s = 1:numel(sheetsToTry)
                try
                    tmp = readmatrix(pathOrFile,'Sheet',sheetsToTry{s});
                    if ~isempty(tmp) && isnumeric(tmp)
                        M = tmp; src = sprintf('%s (sheet=%s)', pathOrFile, sheetsToTry{s}); return;
                    end
                catch, end
            end
            try tmp = readmatrix(pathOrFile); if ~isempty(tmp) && isnumeric(tmp), M=tmp; src=pathOrFile; return; end; end
        elseif strcmpi(ext,'.mat')
            S = load(pathOrFile);
            cand = {'A','A_true','Adj','adj','groundTruth'};
            fn = fieldnames(S);
            for k=1:numel(cand)
                if ismember(cand{k}, fn)
                    M = S.(cand{k}); src = sprintf('%s (var=%s)', pathOrFile, cand{k}); return;
                end
            end
            if numel(fn)==1 && ismatrix(S.(fn{1}))
                M = S.(fn{1}); src = sprintf('%s (var=%s)', pathOrFile, fn{1}); return;
            end
        end
    catch
    end
end

function B = robust_read_pred(pathOrFile, candidateSheets)
    B = [];
    if isempty(pathOrFile), return; end
    [~,~,ext] = fileparts(pathOrFile);
    try
        if strcmpi(ext,'.xlsx') || strcmpi(ext,'.xls')
            for s = 1:numel(candidateSheets)
                try
                    tmp = readmatrix(pathOrFile,'Sheet',candidateSheets{s});
                    if ~isempty(tmp) && isnumeric(tmp)
                        B = tmp; return;
                    end
                catch, end
            end
            try tmp = readmatrix(pathOrFile); if ~isempty(tmp) && isnumeric(tmp), B = tmp; return; end; end
        elseif strcmpi(ext,'.mat')
            S = load(pathOrFile);
            fn = fieldnames(S);
            for k=1:numel(candidateSheets)
                cs = candidateSheets{k};
                if ismember(cs, fn)
                    B = S.(cs); return;
                end
            end
            for i=1:numel(fn)
                v = S.(fn{i});
                if isnumeric(v) && ismatrix(v) && size(v,1)>1 && size(v,2)>1
                    B = v; return;
                end
            end
        end
    catch
    end
end

function Mout = ensure_square_NxN(Min, Ntarget)
    if isempty(Min)
        Mout = []; return;
    end
    if ~isnumeric(Min), error('Input must be numeric.'); end
    [r,c] = size(Min);
    if r==Ntarget && c==Ntarget
        Mout = Min;
    elseif isvector(Min) && numel(Min) >= Ntarget^2
        tmp = Min(:);
        Mout = reshape(tmp(1:Ntarget^2), Ntarget, Ntarget);
    elseif r*c >= Ntarget^2
        tmp = Min(:);
        Mout = reshape(tmp(1:Ntarget^2), Ntarget, Ntarget);
    else
        error('Cannot normalize matrix %dx%d to %dx%d', r,c,Ntarget,Ntarget);
    end
    Mout(~isfinite(Mout)) = 0;
    Mout(1:Ntarget+1:end) = 0;
end

function row = pack_row(Nval,Tval,dtval,sparsityval,snr_dbval,seedval,wval,sigma_rval,tauval,rho_targetval,smooth_sigmaval,out)
    row = [Nval,Tval,dtval,sparsityval,snr_dbval,seedval,wval,sigma_rval,tauval,rho_targetval,smooth_sigmaval, ...
        out.auprc_pos, out.auprc_neg, out.auroc_pos, out.auroc_neg, out.auprc_pos_CI(1), out.auprc_pos_CI(2), ...
        out.auprc_neg_CI(1), out.auprc_neg_CI(2), out.r_pos, out.r_neg, out.r_overall, out.r_pred_on_predicted_edges, ...
        out.mae_pos, out.mae_neg, out.mae_overall, out.rmse_pos, out.rmse_neg, out.rmse_overall, ...
        out.sign_acc_pos, out.sign_acc_neg, out.sign_acc_overall, out.true_density, out.pred_density, out.SHD, out.SHD_signed, ...
        out.precision_at_match, out.recall_at_match, out.F1_at_match];
end

function out = compare_connectivity_matrices(A, B, nboot, tol)
    if nargin<4, tol=1e-6; end
    if nargin<3, nboot=0; end
    N = size(A,1);
    offdiag = ~eye(N);
    vecA = A(offdiag);
    vecB = B(offdiag);
    labels_pos = double(vecA>0);
    labels_neg = double(vecA<0);
    scores = vecB;
    if numel(unique(labels_pos))>1
        [~,~,~,auroc_pos] = perfcurve(labels_pos, scores, 1);
        [~,~,~,auprc_pos] = perfcurve(labels_pos, scores, 1, 'xCrit','reca','yCrit','prec');
    else; auroc_pos=NaN; auprc_pos=NaN; end
    if numel(unique(labels_neg))>1
        [~,~,~,auroc_neg] = perfcurve(labels_neg, -scores, 1);
        [~,~,~,auprc_neg] = perfcurve(labels_neg, -scores, 1, 'xCrit','reca','yCrit','prec');
    else; auroc_neg=NaN; auprc_neg=NaN; end
    mask_pos = (A>0) & offdiag;
    mask_neg = (A<0) & offdiag;
    if any(mask_pos(:))
        r_pos = corr(B(mask_pos), A(mask_pos));
        mae_pos = mean(abs(B(mask_pos) - A(mask_pos)));
        rmse_pos = sqrt(mean((B(mask_pos) - A(mask_pos)).^2));
        sign_acc_pos = mean(sign(B(mask_pos)) == sign(A(mask_pos)));
    else; r_pos=NaN; mae_pos=NaN; rmse_pos=NaN; sign_acc_pos=NaN; end
    if any(mask_neg(:))
        r_neg = corr(B(mask_neg), A(mask_neg));
        mae_neg = mean(abs(B(mask_neg) - A(mask_neg)));
        rmse_neg = sqrt(mean((B(mask_neg) - A(mask_neg)).^2));
        sign_acc_neg = mean(sign(B(mask_neg)) == sign(A(mask_neg)));
    else; r_neg=NaN; mae_neg=NaN; rmse_neg=NaN; sign_acc_neg=NaN; end
    r_overall = corr(vecA, vecB, 'rows','complete');
    mae_overall = mean(abs(vecB-vecA));
    rmse_overall = sqrt(mean((vecB-vecA).^2));
    sign_acc_overall = mean(sign(vecB)==sign(vecA));
    true_density = mean(vecA~=0);
    pred_density = mean(abs(vecB)>tol);
    pred_adj = false(N);
    nz = abs(B)>tol; nz(eye(N)==1)=false;
    if any(nz(:))
        pred_adj(nz)=true;
    else
        Kpos = sum(labels_pos==1); Kneg = sum(labels_neg==1);
        [~, idx_desc] = sort(scores,'descend');
        [~, idx_asc]  = sort(scores,'ascend');
        [rows_off, cols_off] = find(offdiag);
        if Kpos>0
            sel_pos = idx_desc(1:min(Kpos,numel(idx_desc)));
            pred_adj(sub2ind([N,N], rows_off(sel_pos), cols_off(sel_pos))) = true;
        end
        if Kneg>0
            sel_neg = idx_asc(1:min(Kneg,numel(idx_asc)));
            pred_adj(sub2ind([N,N], rows_off(sel_neg), cols_off(sel_neg))) = true;
        end
    end
    true_adj = (A~=0);
    SHD = sum(xor(true_adj(offdiag), pred_adj(offdiag)));
    adj_mismatch = sum(xor(true_adj(offdiag), pred_adj(offdiag)));
    common_mask = true_adj & pred_adj;
    signed_mismatch = 0;
    if any(common_mask(:))
        signed_mismatch = sum(sign(B(common_mask)) ~= sign(A(common_mask)));
    end
    SHD_signed = adj_mismatch + signed_mismatch;
    pred_binary_vec = pred_adj(offdiag);
    true_binary_vec = true_adj(offdiag);
    tp = sum(pred_binary_vec & true_binary_vec);
    fp = sum(pred_binary_vec & ~true_binary_vec);
    fn = sum(~pred_binary_vec & true_binary_vec);
    precision = tp / (tp + fp + eps);
    recall = tp / (tp + fn + eps);
    F1 = 2 * precision * recall / (precision + recall + eps);
    if any(pred_adj(offdiag))
        r_pred_on_predicted_edges = corr(vecA(pred_adj(offdiag)), vecB(pred_adj(offdiag)), 'rows','complete');
    else
        r_pred_on_predicted_edges = NaN;
    end
    if nboot>0
        rng(0);
        n_off = numel(vecA);
        auprc_pos_bs = nan(nboot,1);
        auprc_neg_bs = nan(nboot,1);
        for b=1:nboot
            idx = randsample(n_off, n_off, true);
            lb = labels_pos(idx);
            sb = scores(idx);
            if numel(unique(lb))>1
                [~,~,~,auprc_pos_bs(b)] = perfcurve(lb,sb,1,'xCrit','reca','yCrit','prec');
            else
                auprc_pos_bs(b) = NaN;
            end
            ln = labels_neg(idx);
            if numel(unique(ln))>1
                [~,~,~,auprc_neg_bs(b)] = perfcurve(ln,-sb,1,'xCrit','reca','yCrit','prec');
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
    out.auprc_pos = auprc_pos; out.auprc_neg = auprc_neg;
    out.auroc_pos = auroc_pos; out.auroc_neg = auroc_neg;
    out.auprc_pos_CI = auprc_pos_CI; out.auprc_neg_CI = auprc_neg_CI;
    out.r_pos = r_pos; out.r_neg = r_neg; out.r_overall = r_overall; out.r_pred_on_predicted_edges = r_pred_on_predicted_edges;
    out.mae_pos = mae_pos; out.mae_neg = mae_neg; out.mae_overall = mae_overall;
    out.rmse_pos = rmse_pos; out.rmse_neg = rmse_neg; out.rmse_overall = rmse_overall;
    out.sign_acc_pos = sign_acc_pos; out.sign_acc_neg = sign_acc_neg; out.sign_acc_overall = sign_acc_overall;
    out.true_density = true_density; out.pred_density = pred_density;
    out.SHD = SHD; out.SHD_signed = SHD_signed;
    out.precision_at_match = precision; out.recall_at_match = recall; out.F1_at_match = F1;
end

function v = try_parse_key_number(keys, prefix)
    v = NaN;
    for i=1:numel(keys)
        k = keys{i};
        if startsWith(lower(k), lower(prefix))
            s = regexp(k, [prefix '([0-9\.,]+)'], 'tokens', 'once');
            if ~isempty(s)
                % restore decimal dot if comma used
                sstr = strrep(s{1}, ',', '.');
                v = str2double(sstr);
                return;
            end
        end
    end
end

function wid = getWorkerID()
try
    wid = getCurrentTask().ID;
catch
    wid = 0;
end
end