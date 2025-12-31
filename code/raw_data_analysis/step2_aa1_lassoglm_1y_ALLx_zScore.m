% aa1_lassoglm_1y_ALLx_zScore.m
% Process all groundTruth_* files in ../data/groundTruth.
% Parallelizes over files (one file -> one worker). Per-file GLM loop is serial
% to avoid nested parfor issues.
%
% Optional (set in base before running):
%   seedList = [42 101 ...];   % process only matching seeds
%   sheetType = 'Z';           % sheet to read (default 'Z')
%   distributionSelect = 'normal';

warning('off','all');

%% --- locate input/output folders relative to script ---
fullpath = mfilename('fullpath');
if isempty(fullpath)
    fullpath = matlab.desktop.editor.getActiveFilename; % fallback in editor
end
[scriptFolder, scriptName, ~] = fileparts(fullpath);
parentFolder = fileparts(scriptFolder);

% canonical data/result folders
pathIn   = fullfile(parentFolder, 'data', 'groundTruth');
pathOut  = fullfile(parentFolder, 'data', scriptName);
if ~isfolder(pathIn)
    error('Input folder not found: %s', pathIn);
end
if ~exist(pathOut,'dir'), mkdir(pathOut); end

% GLM defaults (can be defined in base before calling run)
if ~exist('sheetType','var'), sheetType = 'Z'; end
if ~exist('distributionSelect','var'), distributionSelect = 'normal'; end

% optional seed filtering (seedList or single seedval may be present in base)
seedListFilter = [];
if exist('seedList','var') && ~isempty(seedList)
    seedListFilter = double(seedList(:))';
elseif exist('seedval','var') && ~isempty(seedval)
    seedListFilter = double(seedval);
end

%% --- discover files ---
d = dir(fullfile(pathIn,'groundTruth_*.xlsx'));
files = {d.name};
if isempty(files)
    error('No groundTruth_*.xlsx files found in %s', pathIn);
end

% apply seed filter if present
if ~isempty(seedListFilter)
    keep = false(size(files));
    for k=1:numel(files)
        tok = regexp(files{k},'seed[_-]?([0-9]+)','tokens','once');
        if ~isempty(tok) && any(str2double(tok{1}) == seedListFilter)
            keep(k) = true;
        end
    end
    files = files(keep);
    if isempty(files)
        error('No groundTruth files matched seedList filter.');
    end
end

fprintf('Found %d groundTruth files to process.\n', numel(files));

%% --- prepare pool and worker paths ---
pool = gcp('nocreate');
if isempty(pool)
    try
        pool = parpool; %#ok<NASGU>
    catch
        warning('Unable to start parpool automatically. Will run serially.');
        pool = [];
    end
end
if ~isempty(pool)
    pctRunOnAll addpath(scriptFolder);
end

%% --- parallel over files (one file per worker) ---
numFiles = numel(files);
parfor fi = 1:numFiles
    fname = files{fi};
    infile = fullfile(pathIn, fname);

    try
        % set RNG for reproducibility using seed encoded in filename if present
        tok = regexp(fname,'seed[_-]?([0-9]+)','tokens','once');
        if ~isempty(tok)
            sseed = str2double(tok{1});
            if ~isnan(sseed)
                rng(sseed,'twister');
            end
        end
    catch
        % ignore RNG set failures
    end

    try
        fprintf('Worker %d processing %s\n', getWorkerID(), fname);
    catch
        fprintf('Processing %s\n', fname);
    end

    % read Z sheet and transpose to match previous layout
    try
        data = readmatrix(infile,'Sheet',sheetType);
    catch ME
        warning('Failed reading %s sheet %s. Skipping file. (%s)', infile, sheetType, ME.message);
        continue;
    end
    if isempty(data)
        warning('Empty sheet for %s. Skipping.', infile);
        continue;
    end
    dataCa = data';            % columns are neurons, rows timepoints
    [numRows, numCols] = size(dataCa);

    % allocate result matrices (n x n)
    b1_CV1       = NaN(numCols, numCols);
    b0_CV1       = NaN(numCols, numCols);
    lambda_CV1   = NaN(numCols, numCols);

    b1_CV2       = NaN(numCols, numCols);
    b0_CV2       = NaN(numCols, numCols);
    lambda_CV2   = NaN(numCols, numCols);

    mse_CV       = NaN(numCols, numCols);
    Rsqr_CV      = NaN(numCols, numCols);

    b1_AIC       = NaN(numCols, numCols);
    b0_AIC       = NaN(numCols, numCols);

    b1_BIC       = NaN(numCols, numCols);
    b0_BIC       = NaN(numCols, numCols);

    b1_Custom    = NaN(numCols, numCols);
    b0_Custom    = NaN(numCols, numCols);
    mse_Custom   = NaN(numCols, numCols);
    Rsqr_Custom  = NaN(numCols, numCols);

    % iterate neurons (serial inside worker to avoid nested parfor)
    for neuron1 = 1:numCols
        y = dataCa(:, neuron1);
        predIdx = setdiff(1:numCols, neuron1);
        X = dataCa(:, predIdx);

        % CV path via lassoglm
        try
            [bCV, statsCV] = lassoglm(X, y, distributionSelect, 'CV', 10, 'Standardize', true);
        catch ME
            warning('lassoglm failed for file %s neuron %d: %s', fname, neuron1, ME.message);
            continue;
        end

        idxMinDev = statsCV.IndexMinDeviance;
        idx1SE    = statsCV.Index1SE;

        coefMin = NaN(numCols,1);  coefMin(predIdx) = bCV(:, idxMinDev);
        coef1SE = NaN(numCols,1);  coef1SE(predIdx) = bCV(:, idx1SE);

        interceptCV_1 = repmat(statsCV.Intercept(idxMinDev), numCols, 1); interceptCV_1(neuron1) = NaN;
        interceptCV_2 = repmat(statsCV.Intercept(idx1SE), numCols, 1);    interceptCV_2(neuron1) = NaN;

        bestLambdaDev_col = repmat(statsCV.LambdaMinDeviance, numCols, 1); bestLambdaDev_col(neuron1) = NaN;
        bestLambda1SE_col = repmat(statsCV.Lambda1SE, numCols, 1);          bestLambda1SE_col(neuron1) = NaN;

        mse_all = statsCV.Deviance ./ numRows;
        bestMSE_val = mse_all(idxMinDev);
        bestMSE_col = repmat(bestMSE_val, numCols, 1); bestMSE_col(neuron1) = NaN;

        ymean = mean(y);
        TSS = sum((y - ymean).^2);
        RSS = statsCV.Deviance(idxMinDev);
        if TSS ~= 0
            R2_val = 1 - RSS / TSS;
        else
            R2_val = NaN;
        end
        R2_CV_col = repmat(R2_val, numCols, 1); R2_CV_col(neuron1) = NaN;

        % AIC / BIC from CV path
        AIC = statsCV.Deviance + 2 * statsCV.DF;
        BIC = statsCV.Deviance + log(numRows) * statsCV.DF;
        [~, idxAIC] = min(AIC);
        [~, idxBIC] = min(BIC);
        coefAIC = NaN(numCols,1); coefAIC(predIdx) = bCV(:, idxAIC);
        coefBIC = NaN(numCols,1); coefBIC(predIdx) = bCV(:, idxBIC);
        intercept_AIC_col = repmat(statsCV.Intercept(idxAIC), numCols, 1); intercept_AIC_col(neuron1) = NaN;
        intercept_BIC_col = repmat(statsCV.Intercept(idxBIC), numCols, 1); intercept_BIC_col(neuron1) = NaN;

        % Custom holdout
        try
            customValid = cvpartition(size(X,1), 'HoldOut', 0.2);
            Xtrain = X(training(customValid), :);
            ytrain = y(training(customValid));
            Xtest  = X(test(customValid), :);
            ytest  = y(test(customValid));

            [bCustom, statsCustom] = lassoglm(Xtrain, ytrain, distributionSelect, ...
                                              'Lambda', logspace(-4, 1, 100),'Standardize',true);
            interceptRow = reshape(statsCustom.Intercept, 1, []);
            preds = Xtest * bCustom + interceptRow;
            mseCustom = mean((preds - ytest).^2, 1);
            [mseMin, idxCustom] = min(mseCustom);
            coefCustom = NaN(numCols,1); coefCustom(predIdx) = bCustom(:, idxCustom);
            intercept_custom_col = repmat(statsCustom.Intercept(idxCustom), numCols, 1); intercept_custom_col(neuron1) = NaN;
            mseCustomMin_col = repmat(mseMin, numCols, 1); mseCustomMin_col(neuron1) = NaN;

            ymean_test = mean(ytest);
            TSS_test = sum((ytest - ymean_test).^2);
            RSS_test = mseMin * numel(ytest);
            if TSS_test ~= 0
                R2_Custom_val = 1 - RSS_test / TSS_test;
            else
                R2_Custom_val = NaN;
            end
            R2_Custom_col = repmat(R2_Custom_val, numCols, 1); R2_Custom_col(neuron1) = NaN;
        catch
            % in case custom holdout fails, continue with NaNs
            coefCustom = NaN(numCols,1);
            intercept_custom_col = NaN(numCols,1);
            mseCustomMin_col = NaN(numCols,1);
            R2_Custom_col = NaN(numCols,1);
        end

        % assign to global matrices
        b1_CV1(:, neuron1)      = coefMin;
        b0_CV1(:, neuron1)      = interceptCV_1;
        lambda_CV1(:, neuron1)  = bestLambdaDev_col;

        b1_CV2(:, neuron1)      = coef1SE;
        b0_CV2(:, neuron1)      = interceptCV_2;
        lambda_CV2(:, neuron1)  = bestLambda1SE_col;

        mse_CV(:, neuron1)      = bestMSE_col;
        Rsqr_CV(:, neuron1)     = R2_CV_col;

        b1_AIC(:, neuron1)      = coefAIC;
        b0_AIC(:, neuron1)      = intercept_AIC_col;

        b1_BIC(:, neuron1)      = coefBIC;
        b0_BIC(:, neuron1)      = intercept_BIC_col;

        b1_Custom(:, neuron1)   = coefCustom;
        b0_Custom(:, neuron1)   = intercept_custom_col;
        mse_Custom(:, neuron1)  = mseCustomMin_col;
        Rsqr_Custom(:, neuron1) = R2_Custom_col;
    end % neuron loop

    % write per-file outputs
    try
        [~,base,ext] = fileparts(fname);
        outFile = fullfile(pathOut, [scriptName,'_', base, ext]);
        writematrix(b1_CV1,outFile,'Sheet','b1_CV1');
        writematrix(b0_CV1,outFile,'Sheet','b0_CV1');
        writematrix(lambda_CV1,outFile,'Sheet','lambda_CV1');
        writematrix(b1_CV2,outFile,'Sheet','b1_CV2');
        writematrix(b0_CV2,outFile,'Sheet','b0_CV2');
        writematrix(lambda_CV2,outFile,'Sheet','lambda_CV2');
        writematrix(mse_CV,outFile,'Sheet','mse_CV');
        writematrix(Rsqr_CV,outFile,'Sheet','Rsqr_CV');
        writematrix(b1_AIC,outFile,'Sheet','b1_AIC');
        writematrix(b0_AIC,outFile,'Sheet','b0_AIC');
        writematrix(b1_BIC,outFile,'Sheet','b1_BIC');
        writematrix(b0_BIC,outFile,'Sheet','b0_BIC');
        writematrix(b1_Custom,outFile,'Sheet','b1_Custom');
        writematrix(b0_Custom,outFile,'Sheet','b0_Custom');
        writematrix(mse_Custom,outFile,'Sheet','mse_Custom');
        writematrix(Rsqr_Custom,outFile,'Sheet','Rsqr_Custom');
        fprintf('Wrote results: %s\n', outFile);
    catch ME
        warning('Failed writing %s: %s', outFile, ME.message);
    end

end % parfor files

fprintf('aa1_lassoglm_1y_ALLx_zScore: done.\n');

%% helper to get worker id
function wid = getWorkerID()
try
    wid = getCurrentTask().ID;
catch
    wid = 0;
end
end