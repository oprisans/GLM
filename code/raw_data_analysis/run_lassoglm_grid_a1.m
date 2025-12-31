function run_lassoglm_grid_a1(param1Name, param1Vals, param2Name, param2Vals, seeds, outdir)
% RUN_LASSOGLM_GRID  Sweep two parameters and seeds, run full pipeline, aggregate metrics.
% Positional args:
%  run_lassoglm_grid_a1(param1Name, param1Vals, param2Name, param2Vals, seeds, outdir)
% Example:
%  run_lassoglm_grid_a1('sparsityval', 0.05, 'snr_dbval', 0, 1, 'test_out')
%Single seed (fast test):
%run_lassoglm_grid_a1('sparsityval', 0.05:0.05:0.30, 'snr_dbval', 0:5:10, 1, 'test_out')
%Multiple seeds (example 1..5):
%run_lassoglm_grid_a1('sparsityval', 0.05:0.05:0.30, 'snr_dbval', 0:5:10, 1:5, 'test_out')

%% --- resolve project root (use installed code file location) ---
w = which('driver_simulate_var_calcium');
if isempty(w)
    warning(['Cannot find driver_simulate_var_calcium on MATLAB path. ' ...
             'Falling back to current working directory as project root. ' ...
             'If results are not found adjust your path or move to project root.']);
    projectRoot = pwd;
else
    codeDir = fileparts(w);            % .../simulations2/code
    projectRoot = fileparts(codeDir);  % .../simulations2
end
fprintf('Project root resolved to: %s\n', projectRoot);

%% --- outdir handling (make absolute under project root if relative) ---
if nargin<6 || isempty(outdir), outdir = 'batch_out'; end
if ~is_absolute_path(outdir)
    outdir = fullfile(projectRoot, outdir);
end
if ~isfolder(outdir), mkdir(outdir); end
fprintf('Using outdir: %s\n', outdir);

% --- script/function names (adjust if your filenames differ) ---
script_code1 = 'driver_simulate_var_calcium';
script_code2 = 'a1_lassoglm_1y_1x_zScore';
script_code3 = 'compare_groundTruth_vs_a1_lassoglm_1x_1y_zScore';

% --- search patterns and folders (STRICT: only known result folders) ---
ground_truth_pattern = 'groundTruth_*';
lasso_out_pattern     = 'a1_lassoglm_1y_1x_zScore*';
metrics_pattern_glob  = '*comp_groundTruth*';

searchDirs = { fullfile(projectRoot,'data','groundTruth'), ...
               fullfile(projectRoot,'data','a1_lassoglm_1y_1x_zScore'), ...
               fullfile(projectRoot,'results'), ...
               fullfile(projectRoot,'results','zScore') };

allowedExts = {'.xls', '.xlsx', '.mat', '.csv'};  % acceptable metric file types

% baseline fixed params (defaults; adjust if desired)
assignin('base','Nval',10);
assignin('base','Tval',2400);
assignin('base','dtval',0.05);
assignin('base','sparsityval',0.05);
assignin('base','snr_dbval',0);
assignin('base','wval',0.25);
assignin('base','sigma_rval',0.2);
assignin('base','rho_targetval',0.9);
assignin('base','tauval',1.2);
assignin('base','smooth_sigmaval',0.2);

% Prepare lists to record discovered files (we do not move originals)
gtFilesFound = {};
lassoFilesFound = {};
metricsFilesFound = {};

% iterate grid
totalRuns = numel(param1Vals)*numel(param2Vals)*numel(seeds);
runCounter = 0;
for i1 = 1:numel(param1Vals)
    v1 = param1Vals(i1);
    for i2 = 1:numel(param2Vals)
        v2 = param2Vals(i2);
        for s = seeds
            runCounter = runCounter + 1;
            fprintf('\nRun %d/%d : %s=%g, %s=%g, seed=%d\n', runCounter, totalRuns, param1Name, v1, param2Name, v2, s);

            % set parameters in base workspace
            assignin('base', param1Name, v1);
            assignin('base', param2Name, v2);
            assignin('base','seedval', s);
            rng(s);

            % --- run script 1 (simulate) ---
            try run_or_feval(script_code1); catch ME, warning('code1 error: %s', ME.message); end
            pause(0.1);
            gtFile = find_first_matching_file(ground_truth_pattern, searchDirs);
            if ~isempty(gtFile), gtFilesFound{end+1} = gtFile; end

            % --- run script 2 (lasso) ---
            try run_or_feval(script_code2); catch ME, warning('code2 error: %s', ME.message); end
            pause(0.1);
            lassoFile = find_first_matching_file(lasso_out_pattern, searchDirs);
            if ~isempty(lassoFile), lassoFilesFound{end+1} = lassoFile; end

            % --- run script 3 (compare -> metrics) ---
            try run_or_feval(script_code3); catch ME, warning('code3 error: %s', ME.message); end
            pause(0.1);
            metricsFileFound = find_first_matching_file(metrics_pattern_glob, searchDirs);
            if ~isempty(metricsFileFound)
                metricsFilesFound{end+1} = metricsFileFound;
            else
                warning('metrics file not found for this run.');
            end
        end
    end
end

% ----------------- Aggregate metrics in-place (no moving) -----------------
% Prefer canonical zScore results folder for output CSV if present
preferredResultsDir = fullfile(projectRoot,'results','zScore');
if isfolder(preferredResultsDir)
    targetAggregateFolder = preferredResultsDir;
else
    targetAggregateFolder = outdir; % fallback
end

% Discover candidate metric files across configured searchDirs (do not touch code folder)
candidateFiles = {};
for k = 1:numel(searchDirs)
    sd = searchDirs{k};
    if ~isfolder(sd), continue; end
    dd = dir(fullfile(sd, '*comp_groundTruth*.*'));
    for j = 1:numel(dd)
        candidateFiles{end+1,1} = fullfile(dd(j).folder, dd(j).name); %#ok<SAGROW>
    end
end
% include any files recorded during runs
candidateFiles = [candidateFiles; metricsFilesFound(:)];
candidateFiles = unique(candidateFiles);

if isempty(candidateFiles)
    warning('No compare/metrics files found in configured result folders. Nothing to aggregate.');
    return;
end

% Read candidate files into tables
tables = {};
for k = 1:numel(candidateFiles)
    fname = candidateFiles{k};
    if ~isfile(fname)
        warning('Candidate file not found (skipping): %s', fname);
        continue;
    end
    [~,n,ext] = fileparts(fname); ext = lower(ext);
    if ~ismember(ext, allowedExts)
        warning('Skipping file with unsupported extension: %s', fname);
        continue;
    end
    try
        if ismember(ext, {'.xls', '.xlsx', '.csv'})
            T = readtable(fname,'Sheet',1);
        elseif strcmp(ext,'.mat')
            S = load(fname);
            fn = fieldnames(S);
            T = [];
            for fni = 1:numel(fn)
                v = S.(fn{fni});
                if istable(v)
                    T = v; break;
                elseif isstruct(v)
                    try T = struct2table(v); break; catch, end
                end
            end
            if isempty(T)
                warning('No table-like variable in MAT file, skipping: %s', fname);
                continue;
            end
        else
            warning('Unhandled extension for %s; skipping', fname);
            continue;
        end
        if ~istable(T)
            warning('File did not produce a table; skipping: %s', fname);
            continue;
        end
        % Tag the source file base (no path)
        T.source_file = repmat({[n ext]}, height(T), 1);
        tables{end+1} = T; %#ok<SAGROW>
        fprintf('Read: %s (%d rows)\n', fname, height(T));
    catch ME
        warning('Failed reading %s: %s', fname, ME.message);
    end
end

if isempty(tables)
    warning('No readable metric tables found among discovered candidate files.');
    return;
end

% Harmonize columns safely
allNames = {};
for k = 1:numel(tables)
    allNames = [allNames; tables{k}.Properties.VariableNames(:)]; %#ok<AGROW>
end
allNames = unique(allNames);

for k = 1:numel(tables)
    T = tables{k};
    for v = 1:numel(allNames)
        if ~ismember(allNames{v}, T.Properties.VariableNames)
            T.(allNames{v}) = repmat(missing, height(T), 1);
        end
    end
    T = T(:, allNames);
    tables{k} = T;
end

aggT = vertcat(tables{:});
aggCSVpath = fullfile(targetAggregateFolder,'aggregated_metrics.csv');
writetable(aggT, aggCSVpath);
fprintf('Wrote aggregated CSV to: %s (rows=%d)\n', aggCSVpath, height(aggT));

end

%% ----------------- helper functions -----------------
function run_or_feval(name)
% Try to call a function, otherwise run as a script
if exist([name,'.m'],'file')==2 || exist(name,'file')==2
    try
        feval(name);
    catch
        try evalin('base', sprintf('run(''%s'')', name)); catch, error('Failed to run %s', name); end
    end
else
    error('Script/function %s not found on path.', name);
end
end

function fpath = find_first_matching_file(pattern, searchDirs)
% Look for the first match for pattern inside the ordered searchDirs
fpath = '';
for k = 1:numel(searchDirs)
    sd = searchDirs{k};
    if ~isfolder(sd), continue; end
    try
        d = dir(fullfile(sd, pattern));
    catch
        d = [];
    end
    if ~isempty(d)
        fpath = fullfile(sd, d(1).name);
        return;
    end
    % fallback substring match within sd
    dAll = dir(fullfile(sd, '*.*'));
    names = {dAll.name};
    plain = lower(strrep(pattern,'*',''));
    idx = find(contains(lower(names), plain), 1);
    if ~isempty(idx)
        fpath = fullfile(sd, names{idx});
        return;
    end
end
% final fallback: try current folder
d = dir(pattern);
if ~isempty(d)
    fpath = fullfile(pwd, d(1).name);
end
end

function tf = is_absolute_path(p)
% crude absolute-path test for unix/windows
if isempty(p), tf = false; return; end
if ispc
    tf = ~isempty(regexp(p,'^[A-Za-z]:[\\/]', 'once'));
else
    tf = startsWith(p, '/');
end
end