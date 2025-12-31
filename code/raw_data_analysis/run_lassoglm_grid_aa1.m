function run_lassoglm_grid_aa1(param1Name, param1Vals, param2Name, param2Vals, seeds, outdir)
% RUN_LASSOGLM_GRID  Sweep two parameters and seeds, run full pipeline in parallel, aggregate metrics.
% Usage:
%   run_lassoglm_grid('sparsityval', 0.05:0.05:0.30, 'snr_dbval', 0:5:10, 1:10, 'test_out')
%
% Notes:
%  - This version parallelizes over seeds using parfor.
%  - Scripts are executed as scripts (run) inside each worker's base workspace so
%    they see assignin('base',...) variables.
%  - Ensure your code files are on path. Start a parpool prior to calling for faster startup.
%
%
%parpool('local')   % optional: choose #workers
%run_lassoglm_grid_aa1('sparsityval', 0.05:0.05:0.50, 'snr_dbval', 0:5:10, 1:30, 'test_out')

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
script_code2 = 'aa1_lassoglm_1y_ALLx_zScore';
script_code3 = 'compare_groundTruth_vs_aa1_lassoglm_1x_ALLy_zScore';

% Get full file paths for the scripts (so workers can run them directly)
f1 = which(script_code1);
f2 = which(script_code2);
f3 = which(script_code3);
if isempty(f1) || isempty(f2) || isempty(f3)
    warning('One or more script files not found on path. Check that %s, %s, %s are accessible.', script_code1, script_code2, script_code3);
end
% escape single quotes for safe insertion into run('...') string
if ~isempty(f1), f1q = strrep(f1,'''',''''''); else f1q = ''; end
if ~isempty(f2), f2q = strrep(f2,'''',''''''); else f2q = ''; end
if ~isempty(f3), f3q = strrep(f3,'''',''''''); else f3q = ''; end

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

% ensure a parallel pool exists
poolobj = gcp('nocreate');
if isempty(poolobj)
    try
        parpool; % will use default profile and worker count
    catch
        warning('Could not start parpool automatically. Parallel execution will fail if no pool present.');
    end
end

for i1 = 1:numel(param1Vals)
    v1 = param1Vals(i1);
    for i2 = 1:numel(param2Vals)
        v2 = param2Vals(i2);

        % Broadcast v1/v2 into base (main client) as well (optional)
        assignin('base', param1Name, v1);
        assignin('base', param2Name, v2);

        % Run all seeds in parallel. parfor index is integer.
        seedList = seeds(:); % column vector
        nSeeds = numel(seedList);
        fprintf('\nStarting parallel seeds for %s=%g, %s=%g (nSeeds=%d)\n', param1Name, v1, param2Name, v2, nSeeds);

        % Prepare local copies for parfor (broadcast)
        f1_local = f1q;
        f2_local = f2q;
        f3_local = f3q;
        p1 = param1Name;
        p2 = param2Name;

        parfor si = 1:nSeeds
            s = seedList(si);
            % update run counter display is not recommended inside parfor; print minimal info
            fprintf('  worker %d running seed %d (v1=%g v2=%g)\n', getWorkerID(), s, v1, v2);

            % set variables in worker base workspace so scripts can use them
            try
                assignin('base', p1, v1);
                assignin('base', p2, v2);
                assignin('base', 'seedval', s);
                rng(s);
            catch ME
                warning('assignin/rng failed on worker for seed %d: %s', s, ME.message);
            end

            % run scripts in worker base workspace (use full path)
            if ~isempty(f1_local)
                try evalin('base', sprintf('run(''%s'')', f1_local)); catch ME, warning('simulate error seed %d: %s', s, ME.message); end
            else
                warning('simulate script path empty; skipping simulate for seed %d', s);
            end

            if ~isempty(f2_local)
                try evalin('base', sprintf('run(''%s'')', f2_local)); catch ME, warning('lasso error seed %d: %s', s, ME.message); end
            else
                warning('lasso script path empty; skipping lasso for seed %d', s);
            end

            if ~isempty(f3_local)
                try evalin('base', sprintf('run(''%s'')', f3_local)); catch ME, warning('compare error seed %d: %s', s, ME.message); end
            else
                warning('compare script path empty; skipping compare for seed %d', s);
            end
            % don't attempt to collect file names inside parfor — workers write results into results folder
        end % parfor seeds

        % update runCounter for the completed block of seeds
        runCounter = runCounter + nSeeds;
        fprintf('Completed seeds for %s=%g, %s=%g  (totalRunsCompleted ~ %d/%d)\n', param1Name, v1, param2Name, v2, min(runCounter,totalRuns), totalRuns);

        % small pause to allow file system to flush
        pause(0.1);
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

% ---------------- CORRECTED: harmonize columns and combine ----------------
% build list of variable-name cell arrays
nameCells = cellfun(@(tt) tt.Properties.VariableNames(:), tables, 'UniformOutput', false);
allNames = unique(vertcat(nameCells{:}));

% pad missing columns in each table and reorder columns consistently
for k = 1:numel(tables)
    Ttmp = tables{k};
    for v = 1:numel(allNames)
        if ~ismember(allNames{v}, Ttmp.Properties.VariableNames)
            Ttmp.(allNames{v}) = repmat(missing, height(Ttmp), 1);
        end
    end
    Ttmp = Ttmp(:, allNames);
    tables{k} = Ttmp;
end

aggT = vertcat(tables{:});
aggCSVpath = fullfile(targetAggregateFolder,'aggregated_metrics.csv');
writetable(aggT, aggCSVpath);
fprintf('Wrote aggregated CSV to: %s (rows=%d)\n', aggCSVpath, height(aggT));

end

%% ----------------- helper functions -----------------
function wid = getWorkerID()
% returns worker ID if inside parfor, otherwise 0 for client
try
    wid = getCurrentTask().ID;
catch
    wid = 0;
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