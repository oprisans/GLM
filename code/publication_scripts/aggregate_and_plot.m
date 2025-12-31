function aggregate_and_plot(outdir, param1Name, param2Name, metricNames, Nlist, saveFigs)
% AGGREGATE_AND_PLOT  Read aggregated CSV or metrics files, compute mean+95%CI
% and plot metrics vs param1 with one curve per param2 value and subplot per N.
% Usage:
%  aggregate_and_plot('batch_out','sparsityval','snr_dbval', {'auprc_pos','F1_at_match'}, [10 20 30], true)
if nargin<6, saveFigs = true; end

% --- Resolve project root by locating a known code file (robust to cwd) ---
w = which('driver_simulate_var_calcium');
if ~isempty(w)
    codeDir = fileparts(w);
    projectRoot = fileparts(codeDir);
else
    warning('driver_simulate_var_calcium not found on path. Using cwd as project root.');
    projectRoot = pwd;
end

% prefer canonical results/zScore folder
resultsZ = fullfile(projectRoot,'results','zScore');
if isfolder(resultsZ)
    csvFile = fullfile(resultsZ,'aggregated_metrics.csv');
else
    % fallback to provided outdir (make absolute if relative)
    if ~is_absolute_path(outdir), outdir = fullfile(projectRoot,outdir); end
    csvFile = fullfile(outdir,'aggregated_metrics.csv');
end

% --- load or build aggregated table ---
if exist(csvFile,'file')==2
    T = readtable(csvFile);
else
    if isfolder(resultsZ)
        T = run_local_aggregation(resultsZ);
        try writetable(T, fullfile(resultsZ,'aggregated_metrics.csv')); catch, warning('Failed to write aggregated CSV to resultsZ'); end
    else
        if ~isfolder(outdir), error('No aggregated CSV and no results folder to scan.'); end
        T = run_local_aggregation(outdir);
        try writetable(T, csvFile); catch, warning('Failed to write aggregated CSV to outdir'); end
    end
end

% Normalize column names and make them valid identifiers
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% Resolve param columns tolerant to variants like 'sparsityval' vs 'sparsity'
p1 = resolve_param_col(T, param1Name);
p2 = resolve_param_col(T, param2Name);

% ensure metrics exist
nMetric = numel(metricNames);
nboot = 1000;

p2vals = unique(T.(p2));
p1vals = unique(T.(p1));
colors = lines(max(1,numel(p2vals)));

for m = 1:nMetric
    metric = metricNames{m};
    mv = matlab.lang.makeValidName(metric);
    if ~ismember(mv, T.Properties.VariableNames)
        warning('Metric %s not found in table. Skipping.', metric);
        continue;
    end

    figure('Name',sprintf('Metric_%s',metric),'NumberTitle','off','Position',[100 100 1200 600]);
    ncols = numel(Nlist);
    for ni = 1:ncols
        Nval = Nlist(ni);
        subplot(1,ncols,ni); hold on; grid on;
        title(sprintf('%s vs %s (N=%d)', metric, param1Name, Nval));
        xlabel(param1Name); ylabel(metric);

        legendEntries = cell(numel(p2vals),1);
        legendIdx = 0;
        for j = 1:numel(p2vals)
            v2 = p2vals(j);
            % select rows for this N and this p2 value
            if ismember('N', T.Properties.VariableNames)
                mask = (T.(p2) == v2) & (T.N == Nval);
            else
                mask = (T.(p2) == v2);
            end
            if ~any(mask)
                continue;
            end

            % compute mean and bootstrap CI per p1 value
            mu = nan(numel(p1vals),1); ci_lo = nan(numel(p1vals),1); ci_hi = nan(numel(p1vals),1); counts = zeros(numel(p1vals),1);
            for k = 1:numel(p1vals)
                v1 = p1vals(k);
                sel = mask & (T.(p1) == v1);
                vals = T{sel, mv};
                vals = vals(~(ismissing(vals) | isnan(vals)));
                counts(k) = numel(vals);
                if counts(k) == 0
                    mu(k) = NaN; ci_lo(k)=NaN; ci_hi(k)=NaN;
                elseif counts(k) == 1
                    mu(k) = mean(vals); ci_lo(k)=mu(k); ci_hi(k)=mu(k);
                else
                    mu(k) = mean(vals);
                    b = bootstrp(nboot, @mean, vals);
                    ci = prctile(b, [2.5 97.5]);
                    ci_lo(k) = ci(1); ci_hi(k) = ci(2);
                end
            end

            % plotting
            hasVar = counts >= 2;
            xplot = double(p1vals);
            legendIdx = legendIdx + 1;
            if any(hasVar)
                xv = xplot(hasVar); lo = ci_lo(hasVar); hi = ci_hi(hasVar);
                xx = [xv; flipud(xv)]; yy = [lo; flipud(hi)];
                patch(xx, yy, colors(j,:), 'FaceAlpha', 0.15, 'EdgeColor','none');
                plot(xplot, mu, '-o', 'Color', colors(j,:), 'LineWidth',1.4, 'MarkerFaceColor','w');
            else
                plot(xplot, mu, '--o', 'Color', colors(j,:), 'LineWidth',1.2, 'MarkerFaceColor','w');
            end
            % annotate counts
            for k = 1:numel(xplot)
                if ~isnan(mu(k))
                    text(xplot(k), mu(k) - 0.02*max(1,abs(mu(k))), sprintf('n=%d', counts(k)), 'FontSize',8, 'HorizontalAlignment','center');
                end
            end
            legendEntries{legendIdx} = sprintf('%s=%.3g (navg=%g)', param2Name, v2, mean(counts(counts>0)));
        end
        legendEntries = legendEntries(~cellfun('isempty',legendEntries));
        if ~isempty(legendEntries)
            legend(legendEntries,'Location','best');
        end
    end

    if saveFigs
        % save figure beside aggregated CSV
        [aggFolder,~,~] = fileparts(csvFile);
        if isempty(aggFolder), aggFolder = outdir; end
        fname = fullfile(aggFolder, sprintf('metric_%s_vs_%s.png', metric, param1Name));
        saveas(gcf, fname);
        fprintf('Saved figure %s\n', fname);
    end
end

end

%% ---------------- helper: build aggregated table from folder (prefers .xlsx) --------------
function T = run_local_aggregation(folder)
d = dir(fullfile(folder,'*comp_groundTruth*.*'));
if isempty(d), error('No metrics files found in %s', folder); end

% prefer .xlsx when pair exists
map = containers.Map();
for k = 1:numel(d)
    [~,name,ext] = fileparts(d(k).name);
    key = name;
    ext = lower(ext);
    if isKey(map,key)
        prev = map(key); [~,~,prevE] = fileparts(prev); prevE = lower(prevE);
        if strcmp(prevE,'.xlsx'), continue; end
        if strcmp(ext,'.xlsx'), map(key) = fullfile(d(k).folder,d(k).name); end
    else
        map(key) = fullfile(d(k).folder,d(k).name);
    end
end

keysList = keys(map);
tables = {};
for i=1:numel(keysList)
    fname = map(keysList{i});
    [~,basename,ext] = fileparts(fname); ext = lower(ext);
    try
        if ismember(ext,{'.xls','.xlsx','.csv'})
            t = readtable(fname,'Sheet',1);
        elseif strcmp(ext,'.mat')
            S = load(fname);
            fn = fieldnames(S); t = [];
            for fi = 1:numel(fn)
                v = S.(fn{fi});
                if istable(v), t = v; break;
                elseif isstruct(v), try t = struct2table(v); break; catch, end
                end
            end
            if isempty(t), warning('No table-like var in %s; skipping', fname); continue; end
        else
            warning('Skipping unsupported extension: %s', fname); continue;
        end
    catch ME
        warning('Failed reading %s: %s', fname, ME.message);
        continue;
    end
    t.source_file = repmat({[basename ext]}, height(t), 1);
    tables{end+1} = t; %#ok<SAGROW>
end

if isempty(tables), error('No readable metric tables found in %s', folder); end

% harmonize and combine
allNames = {};
for k = 1:numel(tables), allNames = [allNames; tables{k}.Properties.VariableNames(:)]; end
allNames = unique(allNames);
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
T = vertcat(tables{:});
end

%% ---------------- helper: tolerant param name resolver ----------------
function col = resolve_param_col(Tbl, requestedName)
% Try likely variants and return the valid column name present in Tbl
vars = Tbl.Properties.VariableNames;
base = strrep(requestedName,' ','');
candidates = unique({ base, regexprep(base,'_?val$',''), [regexprep(base,'_?val$','') 'val'], [regexprep(base,'_?val$','') '_val'] });
for i = 1:numel(candidates)
    cand = matlab.lang.makeValidName(candidates{i});
    if ismember(cand, vars)
        col = cand; return;
    end
end
% final attempt: validified original
cand = matlab.lang.makeValidName(requestedName);
if ismember(cand, vars), col = cand; return; end
error('Param column "%s" not found. Tried variants: %s. Available columns: %s', requestedName, strjoin(candidates,', '), strjoin(vars,', '));
end

%% ---------------- small util ----------------
function tf = is_absolute_path(p)
if isempty(p), tf = false; return; end
if ispc
    tf = ~isempty(regexp(p,'^[A-Za-z]:[\\/]', 'once'));
else
    tf = startsWith(p, '/');
end
end