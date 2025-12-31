%% driver_plot_aggregate_metrics_multi_sheet.m
clear all; close all; clc;
warning('off');

% ---------------- I/O and paths ----------------
fullpath = mfilename('fullpath');
if isempty(fullpath)
    fullpath = matlab.desktop.editor.getActiveFilename;
end
[scriptFolder, scriptName, ~] = fileparts(fullpath);
parentFolder = fileparts(scriptFolder);

method='aa1_lassoglm_1y_ALLx';
dataType='zScore';
seedRange='seed_1_50';
fileName=[method,'_',dataType,'_comp_groundTruth','_',seedRange,'.xlsx'];

subfolderIn='results';
subfolder_2=[dataType,'_',seedRange];
dataPath=fullfile(parentFolder,subfolderIn,method,subfolder_2);

subfolderOut='figures';
resultsPath   = fullfile(parentFolder,subfolderOut,method,  subfolder_2);
if ~exist(resultsPath, 'dir'), mkdir(resultsPath); end

xlsxRead =fullfile(dataPath,fileName);

% default opts (single-sheet by default). You may override font/DPI here.
opts = struct( ...
  'fixedParams', struct('T',12000,'dt',0.05,'seed',42), ...
  'nboot', 500, ...
  'saveFigs', true, ...
  'sheet', 'b1_CV1', ...        % <-- default, will be overwritten in loop below
  'fontName', 'Arial', ...      % font for exported figures (editable)
  'fontSize', 10, ...           % base font size
  'pngDPI', 600 ...             % raster DPI for PNG
);

% list sheets you want to process
sheets = {'b1_CV1','b1_CV2','b1_AIC','b1_BIC','b1_Custom'};
metrics = {'auprc_pos','F1_at_match','r_overall','SHD'};

% remove seed constraint so plots aggregate across the whole seed block
if isfield(opts,'fixedParams') && isfield(opts.fixedParams,'seed')
    opts.fixedParams = rmfield(opts.fixedParams,'seed');
end

% loop over sheets, setting opts.sheet each time
for si = 1:numel(sheets)
    opts.sheet = sheets{si};      % set single sheet for this iteration
    fprintf('Processing sheet: %s\n', opts.sheet);
    plot_aggregate_metrics(xlsxRead, resultsPath, metrics, 'sparsity', 'snr_db', opts);
end

% ----------------------- function -----------------------
function plot_aggregate_metrics(xlsxFile, outFolder, metrics, xParam, groupParam, opts)
% PLOT_AGGREGATE_METRICS(xlsxFile, outFolder, metrics, xParam, groupParam, opts)
% Reads aggregate XLSX and saves figures into outFolder. Supports opts.sheet as
% either a single sheet name OR a cell array of sheet names (will iterate).

% fallback outFolder
if nargin < 2 || isempty(outFolder)
    outFolder = fileparts(xlsxFile);
    if isempty(outFolder), outFolder = pwd; end
end

% defaults for missing args
if nargin < 3 || isempty(metrics), metrics = {'auprc_pos','F1_at_match','r_overall','SHD'}; end
if nargin < 4 || isempty(xParam), xParam = 'sparsity'; end
if nargin < 5 || isempty(groupParam), groupParam = 'snr_db'; end

% normalize opts
if nargin < 6 || isempty(opts)
    opts = struct();
end
% tolerate opts passed as 1-element cell containing a struct
if iscell(opts) && numel(opts)==1 && isstruct(opts{1})
    opts = opts{1};
end
if ~isstruct(opts)
    error('opts must be a struct (or omitted).');
end

% set safe defaults inside opts
if ~isfield(opts,'nboot') || isempty(opts.nboot), opts.nboot = 500; end
if ~isfield(opts,'saveFigs') || isempty(opts.saveFigs), opts.saveFigs = true; end
if ~isfield(opts,'fixedParams') || isempty(opts.fixedParams), opts.fixedParams = struct(); end
if ~isfield(opts,'colors') || isempty(opts.colors), opts.colors = lines(12); end
if ~isfield(opts,'fontName') || isempty(opts.fontName), opts.fontName = 'Arial'; end
if ~isfield(opts,'fontSize') || isempty(opts.fontSize), opts.fontSize = 10; end
if ~isfield(opts,'pngDPI') || isempty(opts.pngDPI), opts.pngDPI = 600; end

% --- support opts.sheet being a cell array of multiple sheets -> recurse and return
if isfield(opts,'sheet') && ~isempty(opts.sheet)
    % normalize to cellstr
    if isstring(opts.sheet), opts.sheet = cellstr(opts.sheet); end
    if ischar(opts.sheet), opts.sheet = {opts.sheet}; end
    if iscell(opts.sheet) && numel(opts.sheet) > 1
        sheets_local = opts.sheet;
        for s_i = 1:numel(sheets_local)
            opts_i = opts;
            opts_i.sheet = sheets_local{s_i};   % single sheet for recursive call
            plot_aggregate_metrics(xlsxFile, outFolder, metrics, xParam, groupParam, opts_i);
        end
        return;
    end
end

% single-sheet name (or empty)
singleSheet = '';
if isfield(opts,'sheet') && ~isempty(opts.sheet)
    if iscell(opts.sheet), singleSheet = opts.sheet{1}; else singleSheet = opts.sheet; end
end

% ---------- read table ----------
if ~isfile(xlsxFile)
    error('Excel file not found: %s', xlsxFile);
end

if ~isempty(singleSheet)
    try
        T = readtable(xlsxFile, 'Sheet', singleSheet);
    catch ME
        warning('Failed reading sheet "%s": %s\nFalling back to default readtable().', string(singleSheet), ME.message);
        T = readtable(xlsxFile);
    end
else
    T = readtable(xlsxFile);
end

% sanitize column names
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% normalize requested names
metrics = cellfun(@matlab.lang.makeValidName, metrics, 'UniformOutput', false);
xParam = matlab.lang.makeValidName(xParam);
groupParam = matlab.lang.makeValidName(groupParam);

% presence checks
if ~ismember(xParam, T.Properties.VariableNames)
    error('xParam "%s" not found in table.', xParam);
end
if ~ismember(groupParam, T.Properties.VariableNames)
    error('groupParam "%s" not found in table.', groupParam);
end

% keep only metrics that exist
metrics = metrics(ismember(metrics, T.Properties.VariableNames));
if isempty(metrics), error('No valid metrics to plot.'); end

% ---------- apply fixedParams filter ----------
fixedFields = fieldnames(opts.fixedParams);
mask = true(height(T),1);
for i=1:numel(fixedFields)
    f = matlab.lang.makeValidName(fixedFields{i});
    if ~ismember(f, T.Properties.VariableNames)
        error('Fixed param "%s" not in table.', fixedFields{i});
    end
    v = opts.fixedParams.(fixedFields{i});
    col = T.(f);
    if isnumeric(col)
        mask = mask & (abs(col - v) <= max(1e-12, eps(max(1,abs(v)))));
    else
        mask = mask & strcmp(string(col), string(v));
    end
end
T = T(mask,:);
if isempty(T), error('No rows remain after applying fixedParams.'); end

% ---------- unique values ----------
xvals_raw = unique(T.(xParam));
if iscell(xvals_raw) || isstring(xvals_raw)
    xv_num = str2double(string(xvals_raw));
    if all(~isnan(xv_num))
        xvals = unique(xv_num); x_is_numeric = true;
    else
        xvals = unique(categorical(xvals_raw)); x_is_numeric = false;
    end
else
    xvals = unique(T.(xParam)); x_is_numeric = isnumeric(xvals);
end
groupVals = unique(T.(groupParam));

% ensure outFolder
if ~exist(outFolder,'dir'), mkdir(outFolder); end

% facet by N if present
if ismember('N', T.Properties.VariableNames)
    facetVals = unique(T.N);
else
    facetVals = NaN;
end

% ----------------- per-metric plotting (keeps your logic) -----------------
for mi = 1:numel(metrics)
    metric = metrics{mi};
    ciLowName = [metric '_CI_low'];
    ciHighName = [metric '_CI_high'];
    hasPerRowCI = ismember(ciLowName, T.Properties.VariableNames) && ismember(ciHighName, T.Properties.VariableNames);

    for fidx = 1:numel(facetVals)
        facetVal = facetVals(fidx);

        fig = figure('Visible','off','Units','normalized','Position',[0.1 0.1 0.75 0.6]);
        suptitle_str = sprintf('%s vs %s', strrep(metric,'_',' '), strrep(xParam,'_',' '));
        if ~isnan(facetVal), suptitle_str = sprintf('%s  (N=%d)', suptitle_str, facetVal); end
        if exist('sgtitle','file')==2, sgtitle(suptitle_str); end

        nGroups = numel(groupVals);
        ncols = min(3, nGroups);
        nrows = ceil(nGroups / ncols);

        for g = 1:numel(groupVals)
            gv = groupVals(g);
            ax = subplot(nrows, ncols, g); hold(ax,'on'); grid(ax,'on');

            meanY = nan(numel(xvals),1);
            loY = nan(numel(xvals),1);
            hiY = nan(numel(xvals),1);
            counts = zeros(numel(xvals),1);
            true_density_vals = nan(numel(xvals),1);

            for xi = 1:numel(xvals)
                xv = xvals(xi);
                if x_is_numeric
                    sel = abs(T.(xParam) - double(xv)) <= 1e-12;
                else
                    sel = T.(xParam) == xv;
                end
                if isnumeric(T.(groupParam))
                    sel = sel & (abs(T.(groupParam) - double(gv)) <= 1e-12);
                else
                    sel = sel & (T.(groupParam) == gv);
                end
                if ~isnan(facetVal)
                    sel = sel & (T.N == facetVal);
                end
                rows = T(sel,:);
                counts(xi) = height(rows);
                if isempty(rows)
                    meanY(xi)=NaN; loY(xi)=NaN; hiY(xi)=NaN; true_density_vals(xi)=NaN; continue;
                end
                vals = rows.(metric); vals = vals(~(isnan(vals)|ismissing(vals)));
                if isempty(vals)
                    meanY(xi)=NaN; loY(xi)=NaN; hiY(xi)=NaN; true_density_vals(xi)=NaN; continue;
                end
                meanY(xi) = mean(vals);

                if ismember('true_density', rows.Properties.VariableNames)
                    tdcol = rows.true_density; tdcol = tdcol(~(isnan(tdcol)|ismissing(tdcol)));
                    if ~isempty(tdcol), true_density_vals(xi) = mean(tdcol); end
                end

                if hasPerRowCI
                    ci_lo_rows = rows.(ciLowName); ci_hi_rows = rows.(ciHighName);
                    if all(~isnan(ci_lo_rows)) && all(~isnan(ci_hi_rows))
                        mean_lo = mean(ci_lo_rows); mean_hi = mean(ci_hi_rows);
                        loY(xi) = max(0, meanY(xi) - mean_lo);
                        hiY(xi) = max(0, mean_hi - meanY(xi));
                    else
                        if numel(vals) > 1
                            b = bootstrp(opts.nboot, @mean, vals); ci = prctile(b, [2.5 97.5]);
                            loY(xi) = meanY(xi) - ci(1); hiY(xi) = ci(2) - meanY(xi);
                        else
                            loY(xi)=0; hiY(xi)=0;
                        end
                    end
                else
                    if numel(vals) > 1
                        b = bootstrp(opts.nboot, @mean, vals); ci = prctile(b, [2.5 97.5]);
                        loY(xi) = meanY(xi) - ci(1); hiY(xi) = ci(2) - meanY(xi);
                    else
                        loY(xi)=0; hiY(xi)=0;
                    end
                end
            end % xi

            if x_is_numeric, xplot = double(xvals); else xplot = 1:numel(xvals); xticks(ax,xplot); xticklabels(ax,string(xvals)); end
            if all(isnan(meanY)), title(ax, sprintf('%s = %g (no data)', groupParam, gv)); continue; end

            c = opts.colors(mod(g-1,size(opts.colors,1))+1, :);
            valid = ~(isnan(loY) | isnan(hiY) | isnan(meanY));
            if any(valid)
                xv_v = xplot(valid); lo_v = meanY(valid)-loY(valid); hi_v = meanY(valid)+hiY(valid);
                xx = [xv_v; flipud(xv_v)]; yy = [lo_v; flipud(hi_v)];
                patch(ax, xx, yy, c, 'FaceAlpha', 0.15, 'EdgeColor','none');
            end
            plot(ax, xplot, meanY, '-o', 'Color', c, 'LineWidth', 1.5, 'MarkerFaceColor','w');

            % chance line for AUPRC variants using true_density if available
            if contains(lower(metric),'auprc') && any(~isnan(true_density_vals))
                td = true_density_vals;
                % fill NaNs with mean of available TD for plotting continuous line
                if any(~isnan(td))
                    td_filled = td; td_filled(isnan(td_filled)) = mean(td(~isnan(td)));
                    plot(ax, xplot, td_filled, ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
                    % annotate first valid point
                    idxFirst = find(~isnan(td),1,'first');
                    if ~isempty(idxFirst)
                        text(ax, xplot(max(1,idxFirst)), td(idxFirst), sprintf(' chance=%.3g', td(idxFirst)), 'FontSize',8, 'Color',[0.2 0.2 0.2]);
                    end
                end
            end

            title(ax, sprintf('%s = %g (n mean=%.1f)', groupParam, double(gv), mean(counts(counts>0))));
            xlabel(ax, strrep(xParam,'_',' '));
            ylabel(ax, strrep(metric,'_',' '));
            grid(ax,'on');
        end % group loop

        % save per-metric figure
        if opts.saveFigs
            safeMetric = regexprep(metric,'[^A-Za-z0-9]','_');
            if ~isnan(facetVal)
                fnameBase = fullfile(outFolder, sprintf('%s_vs_%s_N_%d', safeMetric, xParam, facetVal));
            else
                fnameBase = fullfile(outFolder, sprintf('%s_vs_%s', safeMetric, xParam));
            end
            sheetSuffix = '';
            if ~isempty(singleSheet), sheetSuffix = ['_' char(matlab.lang.makeValidName(string(singleSheet)))]; end
            % --- use centralized saver that writes SVG, EPS, PNG with font/DPI options ---
            save_figure_outputs(fig, fnameBase, sheetSuffix, opts);
        end
        close(fig);
    end % facet
end % metrics

% ----------------- Combined 2x2 summary figure (AUPRC, F1, SHD, r_pos) -----------------
combinedDesired = {'auprc_pos','F1_at_match','SHD','r_pos'};
combinedMetrics = intersect(combinedDesired, T.Properties.VariableNames, 'stable');
if ~isempty(combinedMetrics)
    for fidx = 1:numel(facetVals)
        facetVal = facetVals(fidx);
        figC = figure('Visible','off','Units','normalized','Position',[0.05 0.05 0.85 0.75]);
        sstr = sprintf('Combined metrics (sheet=%s)', string(singleSheet));
        if ~isnan(facetVal), sstr = sprintf('%s (N=%d)', sstr, facetVal); end
        if exist('sgtitle','file')==2, sgtitle(sstr); end

        for mi = 1:min(4,numel(combinedMetrics))
            metric = combinedMetrics{mi};
            ax = subplot(2,2,mi); hold(ax,'on'); grid(ax,'on');
            title(ax, strrep(metric,'_',' '));
            nGroups = numel(groupVals);
            for g = 1:nGroups
                gv = groupVals(g);
                meanY = nan(numel(xvals),1); loY = nan(numel(xvals),1); hiY = nan(numel(xvals),1); trueDen = nan(numel(xvals),1);
                for xi = 1:numel(xvals)
                    xv = xvals(xi);
                    if x_is_numeric
                        sel = abs(T.(xParam) - double(xv)) <= 1e-12;
                    else
                        sel = T.(xParam) == xv;
                    end
                    if isnumeric(T.(groupParam))
                        sel = sel & (abs(T.(groupParam) - double(gv)) <= 1e-12);
                    else
                        sel = sel & (T.(groupParam) == gv);
                    end
                    if ~isnan(facetVal), sel = sel & (T.N == facetVal); end
                    rows = T(sel,:);
                    if isempty(rows)
                        meanY(xi)=NaN; loY(xi)=NaN; hiY(xi)=NaN; trueDen(xi)=NaN; continue;
                    end
                    vals = rows.(metric); vals = vals(~(isnan(vals)|ismissing(vals)));
                    if isempty(vals), meanY(xi)=NaN; loY(xi)=NaN; hiY(xi)=NaN; trueDen(xi)=NaN; continue; end
                    meanY(xi)=mean(vals);
                    if ismember('true_density', rows.Properties.VariableNames)
                        tdcol = rows.true_density; tdcol = tdcol(~(isnan(tdcol)|ismissing(tdcol)));
                        if ~isempty(tdcol), trueDen(xi)=mean(tdcol); end
                    end
                    % CI
                    ciLo = [metric '_CI_low']; ciHi = [metric '_CI_high'];
                    if ismember(ciLo, rows.Properties.VariableNames) && ismember(ciHi, rows.Properties.VariableNames)
                        lo_ep = mean(rows.(ciLo)); hi_ep = mean(rows.(ciHi));
                        loY(xi) = max(0, meanY(xi) - lo_ep);
                        hiY(xi) = max(0, hi_ep - meanY(xi));
                    else
                        if numel(vals) > 1
                            b = bootstrp(opts.nboot, @mean, vals); ci = prctile(b, [2.5 97.5]);
                            loY(xi) = meanY(xi) - ci(1); hiY(xi) = ci(2) - meanY(xi);
                        else
                            loY(xi)=0; hiY(xi)=0;
                        end
                    end
                end % xi

                if x_is_numeric, xplot = double(xvals); else xplot = 1:numel(xvals); xticks(ax,xplot); xticklabels(ax,string(xvals)); end
                c = opts.colors(mod(g-1,size(opts.colors,1))+1,:);
                valid = ~(isnan(loY)|isnan(hiY)|isnan(meanY));
                if any(valid)
                    xv_v = xplot(valid); lo_v = meanY(valid)-loY(valid); hi_v = meanY(valid)+hiY(valid);
                    xx = [xv_v; flipud(xv_v)]; yy = [lo_v; flipud(hi_v)];
                    patch(ax, xx, yy, c, 'FaceAlpha', 0.12, 'EdgeColor','none');
                end
                plot(ax, xplot, meanY, '-o', 'Color', c, 'LineWidth', 1.4, 'MarkerFaceColor','w');

                % chance line if AUPRC
                if contains(lower(metric),'auprc') && any(~isnan(trueDen))
                    td = trueDen; td(isnan(td)) = mean(td(~isnan(td)));
                    plot(ax, xplot, td, ':', 'Color', [0.2 0.2 0.2], 'LineWidth',1);
                end
            end % groups

            xlabel(ax, strrep(xParam,'_',' ')); ylabel(ax, strrep(metric,'_',' '));
            if numel(groupVals) <= 6, legend(ax, cellstr(string(groupVals)),'Location','best'); end
        end % mi

        % save combined figure
        sheetSafe = '';
        if ~isempty(singleSheet), sheetSafe = ['_' char(matlab.lang.makeValidName(string(singleSheet)))]; end
        if ~isnan(facetVal)
            fnameComb = fullfile(outFolder, sprintf('combined_metrics_N_%d%s', facetVal, sheetSafe));
        else
            fnameComb = fullfile(outFolder, sprintf('combined_metrics%s', sheetSafe));
        end
        save_figure_outputs(figC, fnameComb, '', opts);
        close(figC);
    end
end

fprintf('Done. Figures saved to: %s\n', outFolder);
end

% ----------------- helper: centralized figure exporter -----------------
function save_figure_outputs(figHandle, fnameBase, suffix, opts)
% Saves SVG (editable vector), EPS (vector) and PNG (raster) with sensible defaults.
% - figHandle : figure handle
% - fnameBase : full path base (no extension)
% - suffix    : string appended to base (including leading underscore if desired)
% - opts      : struct with optional fields fontName, fontSize, pngDPI
if nargin < 4, opts = struct(); end
if nargin < 3 || isempty(suffix), suffix = ''; end

% defaults
if ~isfield(opts,'fontName') || isempty(opts.fontName), opts.fontName = 'Arial'; end
if ~isfield(opts,'fontSize') || isempty(opts.fontSize), opts.fontSize = 10; end
if ~isfield(opts,'pngDPI') || isempty(opts.pngDPI), opts.pngDPI = 600; end

svgFile = [fnameBase, suffix, '.svg'];
epsFile = [fnameBase, suffix, '.eps'];
pngFile = [fnameBase, suffix, '.png'];

% 1) enforce consistent fonts and sizes (best-effort)
try
    objs = findall(figHandle, '-property', 'FontName');
    for k = 1:numel(objs)
        try objs(k).FontName = opts.fontName; end
    end
    objs2 = findall(figHandle, '-property', 'FontSize');
    for k = 1:numel(objs2)
        try objs2(k).FontSize = opts.fontSize; end
    end
catch
    % nonfatal
end

% 2) choose an appropriate renderer (painters preferred for vector output)
origRenderer = get(figHandle, 'Renderer');
try
    set(figHandle, 'Renderer', 'painters');
catch
    try set(figHandle, 'Renderer', 'opengl'); catch, end
end
drawnow;

% 3) export SVG (vector) - prefer exportgraphics, fallback to print -dsvg
try
    exportgraphics(figHandle, svgFile, 'ContentType', 'vector');
catch ME
    warning('SVG export via exportgraphics failed: %s. Trying print -dsvg fallback.', ME.message);
    try
        print(figHandle, svgFile, '-dsvg');
    catch ME2
        warning('Fallback SVG print failed: %s', ME2.message);
    end
end

% 4) export EPS (vector). Use exportgraphics then print as compatibility fallback.
try
    exportgraphics(figHandle, epsFile, 'ContentType', 'vector');
    % ensure high-res print as well (helps embedded raster parts)
    try print(figHandle, epsFile, '-depsc2', ['-r' num2str(max(300, opts.pngDPI))]); end
catch ME
    warning('EPS export via exportgraphics failed: %s. Trying print fallback.', ME.message);
    try
        print(figHandle, epsFile, '-depsc2', ['-r' num2str(max(300, opts.pngDPI))]);
    catch ME2
        warning('Fallback EPS print failed: %s', ME2.message);
    end
end

% 5) export PNG (raster) at requested DPI
try
    exportgraphics(figHandle, pngFile, 'Resolution', opts.pngDPI);
catch ME
    warning('PNG export via exportgraphics failed: %s. Trying print fallback.', ME.message);
    try
        print(figHandle, pngFile, ['-r' num2str(opts.pngDPI)], '-dpng');
    catch ME2
        warning('Fallback PNG print failed: %s. Using saveas as last resort.', ME2.message);
        try saveas(figHandle, pngFile); catch, end
    end
end

% restore original renderer if possible
try set(figHandle, 'Renderer', origRenderer); catch, end
end