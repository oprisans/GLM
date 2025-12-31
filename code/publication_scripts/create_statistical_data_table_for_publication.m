% summarize_metrics_for_paper.m
% Single script: driver + generate_numeric_latex_tables + generate_supp_tables
% Edit CONFIG section below and run.

clearvars; close all; clc;
warning('off');

% ================= CONFIG =================
% Paths and filenames (edit)
fullpath = mfilename('fullpath');
if isempty(fullpath)
    fullpath = matlab.desktop.editor.getActiveFilename;
end
[scriptFolder, ~, ~] = fileparts(fullpath);
parentFolder = fileparts(scriptFolder);

method='aa1_lassoglm_1y_ALLx';
dataType='zScore';
seedRange='seed_1_50';
fileName=[method,'_',dataType,'_comp_groundTruth','_',seedRange,'.xlsx'];

subfolderIn='results';
subfolder_2=[method filesep dataType '_' seedRange]; % adjust folder arrangement if needed
dataPath=fullfile(parentFolder, subfolderIn, subfolder_2);

subfolderOut='figures';
resultsPath = fullfile(parentFolder, subfolderOut, subfolder_2);
if ~exist(resultsPath,'dir'), mkdir(resultsPath); end

xlsxFile = fullfile(dataPath, fileName);

% sheets to process (set to {} to auto-detect)
sheets = {'b1_CV1','b1_CV2','b1_AIC','b1_BIC','b1_Custom'}; 

% Which Ns to create numeric tables for (e.g., supplementary)
Ns_to_write = [10 20 30];

% Triplets for the Results table (choose 2-3 representative triplets)
triplets = [ struct('N',30,'sparsity',0.05,'snr_db',0), ...
             struct('N',30,'sparsity',0.20,'snr_db',5), ...
             struct('N',30,'sparsity',0.50,'snr_db',10) ];

% metrics / plotting choices (used by functions)
metrics = {'auprc_pos','F1_at_match','SHD','r_pos'};

% numeric table grid defaults (used by generate_numeric_latex_tables)
sparsities = [0.05 0.20 0.50];
snr_values = [0 5 10];

% bootstrap reps for CI fallback
nboot = 500;

% ============== RUN TASKS ================
% 1) Generate per-sheet numeric LaTeX tables for Ns_to_write (one file per sheet×N)
generate_numeric_latex_tables(xlsxFile, resultsPath, sheets, Ns_to_write, sparsities, snr_values, nboot);

% 2) Generate combined supplementary LaTeX table for specified triplets across sheets
generate_supp_tables(xlsxFile, resultsPath, sheets, triplets, nboot);

% ============== END DRIVER ==============


% ------------------ FUNCTIONS ------------------

function generate_numeric_latex_tables(xlsxFile, outDir, sheets, Ns_to_write, sparsities, snr_values, nboot)
% Write fully-numeric LaTeX tables per sheet × N for the specified grid.
% If sheets is empty or {}, detect sheet names from xlsxFile.

if nargin < 2 || isempty(outDir), outDir = pwd; end
if nargin < 3, sheets = {}; end
if nargin < 4 || isempty(Ns_to_write), Ns_to_write = [10 20]; end
if nargin < 5 || isempty(sparsities), sparsities = [0.05 0.2 0.5]; end
if nargin < 6 || isempty(snr_values), snr_values = [0 5 10]; end
if nargin < 7 || isempty(nboot), nboot = 500; end

if ~exist(outDir,'dir'), mkdir(outDir); end

% get sheet list if needed
if isempty(sheets)
    try
        [~,sheetNames] = xlsfinfo(xlsxFile);
        sheets = sheetNames;
    catch
        error('Unable to determine sheets; provide sheets list in driver.');
    end
end

% common metrics to extract
reqMetrics = {'auprc_pos','F1_at_match','SHD','r_pos'};

for si = 1:numel(sheets)
    sheet = sheets{si};
    try
        T = readtable(xlsxFile, 'Sheet', sheet);
    catch ME
        warning('Skipping sheet "%s": %s', sheet, ME.message);
        continue;
    end
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    for Ni = 1:numel(Ns_to_write)
        Nval = Ns_to_write(Ni);
        rows_struct = []; % accumulate rows for LaTeX

        for isp = 1:numel(sparsities)
            sp = sparsities(isp);
            snr = snr_values(min(isp, numel(snr_values)));

            % selection mask
            mask = true(height(T),1);
            if ismember('N', T.Properties.VariableNames)
                mask = mask & (T.N == Nval);
            else
                mask = false(height(T),1);
            end
            if ismember('sparsity', T.Properties.VariableNames)
                col = T.sparsity;
                if isnumeric(col)
                    mask = mask & (abs(col - sp) <= 1e-8);
                else
                    mask = mask & (strcmp(string(col), string(sp)));
                end
            end
            if ismember('snr_db', T.Properties.VariableNames)
                col = T.snr_db;
                if isnumeric(col)
                    mask = mask & (abs(col - snr) <= 1e-8);
                else
                    mask = mask & (strcmp(string(col), string(snr)));
                end
            end

            sub = T(mask,:);
            r.sheet = sheet; r.N = Nval; r.sparsity = sp; r.snr_db = snr;
            % compute stats for requested metrics
            for mm = 1:numel(reqMetrics)
                mname = matlab.lang.makeValidName(reqMetrics{mm});
                if ~ismember(mname, sub.Properties.VariableNames) || isempty(sub)
                    r.([mname '_mean']) = NaN;
                    r.([mname '_lo'])   = NaN;
                    r.([mname '_hi'])   = NaN;
                    continue;
                end
                vals = sub.(mname);
                vals = vals(~(isnan(vals)|ismissing(vals)));
                % prefer per-row CI endpoints if present
                ciLoName = [mname '_CI_low']; ciHiName = [mname '_CI_high'];
                if ismember(ciLoName, sub.Properties.VariableNames) && ismember(ciHiName, sub.Properties.VariableNames)
                    lo_rows = sub.(ciLoName); hi_rows = sub.(ciHiName);
                    valid_lo = lo_rows(~(isnan(lo_rows)|ismissing(lo_rows)));
                    valid_hi = hi_rows(~(isnan(hi_rows)|ismissing(hi_rows)));
                    if ~isempty(valid_lo) && ~isempty(valid_hi) && ~isempty(vals)
                        r.([mname '_mean']) = mean(vals);
                        r.([mname '_lo']) = mean(valid_lo);
                        r.([mname '_hi']) = mean(valid_hi);
                        continue;
                    end
                end
                % fallback to bootstrap CI over vals
                if isempty(vals)
                    r.([mname '_mean']) = NaN;
                    r.([mname '_lo'])   = NaN;
                    r.([mname '_hi'])   = NaN;
                else
                    r.([mname '_mean']) = mean(vals);
                    if numel(vals) > 1
                        b = bootstrp(nboot, @mean, vals);
                        ci = prctile(b, [2.5 97.5]);
                        r.([mname '_lo']) = ci(1);
                        r.([mname '_hi']) = ci(2);
                    else
                        r.([mname '_lo']) = r.([mname '_mean']);
                        r.([mname '_hi']) = r.([mname '_mean']);
                    end
                end
            end % metrics
            rows_struct = [rows_struct; r]; %#ok<AGROW>
        end % sparsities

        % write LaTeX numeric table for this sheet × N
        safeSheet = regexprep(sheet,'[^A-Za-z0-9]','_');
        texName = fullfile(outDir, sprintf('figure_summary_table_%s_N%d.tex', safeSheet, Nval));
        fid = fopen(texName,'w');
        if fid <= 0
            warning('Could not open %s for writing.', texName);
            continue;
        end
        fprintf(fid, '%% Numeric LaTeX table for sheet=%s N=%d\n', sheet, Nval);
        fprintf(fid, '\\begin{table}[ht]\n\\centering\n');
        fprintf(fid, '\\caption{Summary statistics for %s (N=%d). Mean ± 95%% CI.}\n', sheet, Nval);
        fprintf(fid, '\\begin{tabular}{c c c c c c c c}\n\\hline\n');
        fprintf(fid, 'sheet & N & sparsity & snr\\_db & AUPRC & F1 & SHD & r\\_pos \\\\\n\\hline\n');

        for k = 1:numel(rows_struct)
            R = rows_struct(k);
            % format strings with consistent decimals
            au = format_stat(R.auprc_pos_mean, R.auprc_pos_lo, R.auprc_pos_hi, '%.3f');
            f1 = format_stat(R.F1_at_match_mean, R.F1_at_match_lo, R.F1_at_match_hi, '%.3f');
            shd = format_stat(R.SHD_mean, R.SHD_lo, R.SHD_hi, '%.0f');
            rp  = format_stat(R.r_pos_mean, R.r_pos_lo, R.r_pos_hi, '%.3f');
            fprintf(fid, '%s & %d & %.3g & %d & %s & %s & %s & %s \\\\\n', ...
                sheet, R.N, R.sparsity, R.snr_db, au, f1, shd, rp);
        end
        fprintf(fid, '\\hline\n\\end{tabular}\n\\label{tab:summary_%s_N%d}\n\\end{table}\n', safeSheet, Nval);
        fclose(fid);
        fprintf('Wrote numeric LaTeX table: %s\n', texName);

        % print representative lines to console (3 lines chosen as example; here all rows)
        fprintf('\nRepresentative lines for sheet=%s N=%d (paste verbatim):\n', sheet, Nval);
        for k=1:min(numel(rows_struct), 10)
            R = rows_struct(k);
            au = format_stat(R.auprc_pos_mean, R.auprc_pos_lo, R.auprc_pos_hi, '%.3f');
            f1 = format_stat(R.F1_at_match_mean, R.F1_at_match_lo, R.F1_at_match_hi, '%.3f');
            shd = format_stat(R.SHD_mean, R.SHD_lo, R.SHD_hi, '%.0f');
            rp  = format_stat(R.r_pos_mean, R.r_pos_lo, R.r_pos_hi, '%.3f');
            fprintf('- %s | N=%d sparsity=%.2g snr_db=%d: AUPRC = %s; F1 = %s; SHD = %s; r_pos = %s\n', ...
                sheet, R.N, R.sparsity, R.snr_db, au, f1, shd, rp);
        end
    end % N loop
end % sheets loop

    function s = format_stat(mean_v, lo_v, hi_v, fmt)
        if isempty(mean_v) || isnan(mean_v)
            s = 'NA';
        else
            s = sprintf([fmt ' (' fmt '--' fmt ')'], mean_v, lo_v, hi_v);
        end
    end
end


function generate_supp_tables(xlsxFile, outDir, sheets, triplets, nboot)
% Generates a single supplementary LaTeX table (fully numeric) containing the
% requested triplets for each sheet. Writes out 'supplementary_table.tex' to outDir.

if nargin < 2 || isempty(outDir), outDir = pwd; end
if nargin < 3, sheets = {}; end
if nargin < 4 || isempty(triplets)
    error('triplets must be provided (array of structs with fields N, sparsity, snr_db).');
end
if nargin < 5 || isempty(nboot), nboot = 500; end

if isempty(sheets)
    try
        [~,sheetNames] = xlsfinfo(xlsxFile);
        sheets = sheetNames;
    catch
        error('Provide sheets list or ensure xlsfinfo works on your platform.');
    end
end

if ~exist(outDir,'dir'), mkdir(outDir); end

texName = fullfile(outDir, 'supplementary_table.tex');
fid = fopen(texName,'w');
if fid <= 0, error('Could not open %s for writing.', texName); end

fprintf(fid, '%% Supplementary numeric table: triplets across sheets\n');
fprintf(fid, '\\begin{table}[ht]\\centering\n');
fprintf(fid, '\\caption{Summary (mean ± 95\\%% CI) for selected simulation conditions.}\\begin{tabular}{l c c c c c c c}\\hline\n');
fprintf(fid, 'sheet & N & sparsity & snr\\_db & AUPRC & F1 & SHD & r\\_pos \\\\\\hline\n');

for si = 1:numel(sheets)
    sheet = sheets{si};
    try
        T = readtable(xlsxFile, 'Sheet', sheet);
    catch
        warning('Skipping sheet "%s" (read failed).', sheet);
        continue;
    end
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    for ti = 1:numel(triplets)
        q = triplets(ti);
        mask = true(height(T),1);
        if ismember('N', T.Properties.VariableNames)
            mask = mask & (T.N == q.N);
        else
            mask = false(height(T),1);
        end
        if ismember('sparsity', T.Properties.VariableNames)
            col = T.sparsity;
            if isnumeric(col)
                mask = mask & (abs(col - q.sparsity) <= 1e-8);
            else
                mask = mask & (strcmp(string(col), string(q.sparsity)));
            end
        end
        if ismember('snr_db', T.Properties.VariableNames)
            col = T.snr_db;
            if isnumeric(col)
                mask = mask & (abs(col - q.snr_db) <= 1e-8);
            else
                mask = mask & (strcmp(string(col), string(q.snr_db)));
            end
        end

        rows = T(mask,:);
        if isempty(rows)
            fprintf(fid, '%s & %d & %.3g & %d & -- & -- & -- & -- \\\\\n', sheet, q.N, q.sparsity, q.snr_db);
            fprintf('Note: no data for %s N=%d s=%.3g snr=%d\n', sheet, q.N, q.sparsity, q.snr_db);
            continue;
        end

        % compute stats for metrics
        [au_m, au_lo, au_hi] = compute_metric_rows(rows, 'auprc_pos', nboot);
        [f1_m, f1_lo, f1_hi] = compute_metric_rows(rows, 'F1_at_match', nboot);
        [shd_m, shd_lo, shd_hi] = compute_metric_rows(rows, 'SHD', nboot);
        [rp_m, rp_lo, rp_hi] = compute_metric_rows(rows, 'r_pos', nboot);

        au_s = sprintf('%.3f (%.3f--%.3f)', au_m, au_lo, au_hi);
        f1_s = sprintf('%.3f (%.3f--%.3f)', f1_m, f1_lo, f1_hi);
        shd_s = sprintf('%.0f (%.1f--%.1f)', shd_m, shd_lo, shd_hi);
        rp_s = sprintf('%.3f (%.3f--%.3f)', rp_m, rp_lo, rp_hi);

        fprintf(fid, '%s & %d & %.3g & %d & %s & %s & %s & %s \\\\\n', sheet, q.N, q.sparsity, q.snr_db, au_s, f1_s, shd_s, rp_s);
        % also print a representative console line
        fprintf('- %s | N=%d sparsity=%.3g snr_db=%d: AUPRC = %s; F1 = %s; SHD = %s; r_pos = %s\n', sheet, q.N, q.sparsity, q.snr_db, au_s, f1_s, shd_s, rp_s);
    end
    fprintf(fid, '\\hline\n');
end

fprintf(fid, '\\end{tabular}\\label{tab:supp_summary}\\end{table}\n');
fclose(fid);
fprintf('Wrote supplementary table: %s\n', texName);

    function [m, lo, hi] = compute_metric_rows(tab, metricName, nboot_local)
        m = NaN; lo = NaN; hi = NaN;
        mname = matlab.lang.makeValidName(metricName);
        if ~ismember(mname, tab.Properties.VariableNames)
            return;
        end
        vals = tab.(mname);
        vals = vals(~(isnan(vals)|ismissing(vals)));
        if isempty(vals), return; end
        % per-row CI endpoints preferred
        ciLoName = [mname '_CI_low']; ciHiName = [mname '_CI_high'];
        if ismember(ciLoName, tab.Properties.VariableNames) && ismember(ciHiName, tab.Properties.VariableNames)
            lo_rows = tab.(ciLoName); hi_rows = tab.(ciHiName);
            valid_lo = lo_rows(~(isnan(lo_rows)|ismissing(lo_rows)));
            valid_hi = hi_rows(~(isnan(hi_rows)|ismissing(hi_rows)));
            if ~isempty(valid_lo) && ~isempty(valid_hi)
                m = mean(vals);
                lo = mean(valid_lo);
                hi = mean(valid_hi);
                return;
            end
        end
        m = mean(vals);
        if numel(vals) > 1
            b = bootstrp(nboot_local, @mean, vals);
            ci = prctile(b, [2.5 97.5]);
            lo = ci(1); hi = ci(2);
        else
            lo = m; hi = m;
        end
    end
end