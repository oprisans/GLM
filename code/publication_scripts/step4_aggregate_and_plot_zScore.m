function step4_aggregate_and_plot_zScore(metricNames, Nlist, saveFigs)
% AGGREGATE_AND_PLOT_ZSCORE  Build aggregated_metrics.csv from simulations2/results/zScore
% and plot metrics vs sparsity grouped by snr_db and N.

if nargin<1 || isempty(metricNames), metricNames = {'auprc_pos','F1_at_match'}; end
if nargin<2 || isempty(Nlist), Nlist = unique([10]); end
if nargin<3, saveFigs = true; end

% resolve zScore folder from known code file
w = which('driver_simulate_var_calcium');
if isempty(w)
    projectRoot = pwd;
else
    projectRoot = fileparts(fileparts(w));
end
zdir = fullfile(projectRoot,'results','zScore');
if ~isfolder(zdir), error('zScore folder not found: %s', zdir); end

% collect files, prefer .xlsx when pair exists
d = dir(fullfile(zdir,'*comp_groundTruth*.*'));
if isempty(d), error('No *comp_groundTruth* files found in %s', zdir); end
map = containers.Map();
for k=1:numel(d)
    [~,name,ext] = fileparts(d(k).name); ext = lower(ext);
    key = name;
    if isKey(map,key)
        prev = map(key);
        [~,~,prevE] = fileparts(prev); prevE = lower(prevE);
        % prefer .xlsx if available
        if strcmp(prevE,'.xlsx'), continue; end
        if strcmp(ext,'.xlsx'), map(key) = fullfile(d(k).folder,d(k).name); end
    else
        map(key) = fullfile(d(k).folder,d(k).name);
    end
end

% read selected files into tables
keysList = keys(map);
tables = cell(0,1);
for i=1:numel(keysList)
    fname = map(keysList{i});
    [~,basename,ext] = fileparts(fname); ext = lower(ext);
    try
        if ismember(ext,{'.xls','.xlsx'}) 
            t = readtable(fname,'Sheet',1);
        elseif strcmp(ext,'.csv')
            t = readtable(fname);                   % CSV -> no Sheet
        elseif strcmp(ext,'.mat')
            S = load(fname);
            fn = fieldnames(S); t = [];
            for fi=1:numel(fn)
                v = S.(fn{fi});
                if istable(v), t=v; break;
                elseif isstruct(v), try t = struct2table(v); break; catch, end
                end
            end
            if isempty(t), warning('No table-like var in %s; skipping', fname); continue; end
        else
            warning('Unsupported ext %s for %s; skipping', ext, fname); continue;
        end
    catch ME
        warning('Failed reading %s: %s', fname, ME.message); continue;
    end
    t.source_file = repmat({[basename ext]}, height(t), 1);
    tables{end+1,1} = t;
    fprintf('Read: %s (%d rows)\n', fname, height(t));
end
if isempty(tables), error('No readable metric tables found in %s', zdir); end

% ---------------- harmonize columns and combine ----------------
% collect VariableNames safely
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

% now vertically concatenate safely
aggT = vertcat(tables{:});

% write aggregated CSV into zdir
aggCSV = fullfile(zdir,'aggregated_metrics.csv');
writetable(aggT, aggCSV);
fprintf('Wrote aggregated CSV: %s (rows=%d)\n', aggCSV, height(aggT));

% normalize column names BEFORE resolving params
aggT.Properties.VariableNames = matlab.lang.makeValidName(aggT.Properties.VariableNames);

% tolerant param names: prefer 'sparsity' and 'snr_db'
p1 = resolve_param_col(aggT, 'sparsity');
p2 = resolve_param_col(aggT, 'snr_db');

% plotting defaults
colors = lines(16);
nboot = 1000;

for m = 1:numel(metricNames)
    metric = metricNames{m};
    mv = matlab.lang.makeValidName(metric);
    if ~ismember(mv, aggT.Properties.VariableNames)
        warning('Metric %s not found in aggregated table. Skipping.', metric); continue;
    end

    p2vals = unique(aggT.(p2));
    p1vals = unique(aggT.(p1));
    figure('Name',sprintf('Metric_%s',metric),'NumberTitle','off','Position',[100 100 1200 600]);
    ncols = numel(Nlist);
    for ni=1:ncols
        Nval = Nlist(ni);
        subplot(1,ncols,ni); hold on; grid on;
        title(sprintf('%s vs %s (N=%d)', metric, p1, Nval));
        xlabel(p1); ylabel(metric);
        legendEntries = {};
        li = 0;
        for j=1:numel(p2vals)
            v2 = p2vals(j);
            if ismember('N', aggT.Properties.VariableNames)
                mask = (aggT.(p2) == v2) & (aggT.N == Nval);
            else
                mask = (aggT.(p2) == v2);
            end
            if ~any(mask), continue; end
            mu = nan(numel(p1vals),1); ci_lo = mu; ci_hi = mu; counts = zeros(numel(p1vals),1);
            for k=1:numel(p1vals)
                v1 = p1vals(k);
                sel = mask & (aggT.(p1) == v1);
                vals = aggT{sel, mv};
                vals = vals(~(ismissing(vals) | isnan(vals)));
                counts(k) = numel(vals);
                if counts(k) == 0
                    mu(k)=NaN; ci_lo(k)=NaN; ci_hi(k)=NaN;
                elseif counts(k) == 1
                    mu(k)=mean(vals); ci_lo(k)=mu(k); ci_hi(k)=mu(k);
                else
                    mu(k)=mean(vals);
                    b = bootstrp(nboot, @mean, vals);
                    ci = prctile(b, [2.5 97.5]);
                    ci_lo(k)=ci(1); ci_hi(k)=ci(2);
                end
            end
            li = li + 1;
            xplot = double(p1vals);
            hasVar = counts >= 2;
            if any(hasVar)
                xv = xplot(hasVar); lo = ci_lo(hasVar); hi = ci_hi(hasVar);
                xx = [xv; flipud(xv)]; yy = [lo; flipud(hi)];
                patch(xx, yy, colors(mod(j-1,size(colors,1))+1,:), 'FaceAlpha', 0.15, 'EdgeColor','none');
                plot(xplot, mu, '-o', 'Color', colors(mod(j-1,size(colors,1))+1,:), 'LineWidth',1.4, 'MarkerFaceColor','w');
            else
                plot(xplot, mu, '--o', 'Color', colors(mod(j-1,size(colors,1))+1,:), 'LineWidth',1.2, 'MarkerFaceColor','w');
            end
            for k=1:numel(xplot)
                if ~isnan(mu(k)), text(xplot(k), mu(k) - 0.02*max(1,abs(mu(k))), sprintf('n=%d', counts(k)), 'FontSize',8, 'HorizontalAlignment','center'); end
            end
            legendEntries{li} = sprintf('%s=%.3g (navg=%.1f)', p2, v2, mean(counts(counts>0)));
        end
        if ~isempty(legendEntries), legend(legendEntries,'Location','best'); end
    end

    if saveFigs
        fname = fullfile(zdir, sprintf('metric_%s_vs_%s.png', metric, p1));
        saveas(gcf, fname);
        fprintf('Saved figure: %s\n', fname);
    end
end

end

%% ---------------- helpers ----------------
function col = resolve_param_col(Tbl, requestedName)
vars = Tbl.Properties.VariableNames;
base = strrep(requestedName,' ','');
candidates = unique({ base, regexprep(base,'_?val$',''), [regexprep(base,'_?val$','') 'val'], [regexprep(base,'_?val$','') '_val'] });
for i = 1:numel(candidates)
    cand = matlab.lang.makeValidName(candidates{i});
    if ismember(cand, vars), col = cand; return; end
end
% fallback: if exact not found, return first close match by substring
for i=1:numel(vars)
    if contains(lower(vars{i}), lower(regexprep(base,'_?val$',''))), col = vars{i}; return; end
end
error('Param column "%s" not found. Available: %s', requestedName, strjoin(vars,', '));
end