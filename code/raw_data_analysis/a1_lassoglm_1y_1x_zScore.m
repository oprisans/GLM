%% corrected_univariate_glm_with_AIC_BIC_Custom.m
clear all; close all; clc;
%warning('off');

% parallel pool safe
p = gcp('nocreate');
if isempty(p)
    parpool;
end
pctRunOnAll addpath(genpath(pwd));

% ---------------- I/O and paths ----------------
fullpath = mfilename('fullpath');
if isempty(fullpath)
    fullpath = matlab.desktop.editor.getActiveFilename;
end
[scriptFolder, scriptName, ~] = fileparts(fullpath);
parentFolder = fileparts(scriptFolder);

dataPath   = fullfile(parentFolder, 'data');
resultsPath= fullfile(parentFolder, 'results');

% Input folder containing ground-truth / data files (adjust if needed)
dataType = 'groundTruth';
pathIn   = fullfile(dataPath, dataType);

% GLM settings
GLM_method = 'lassoglm';
sheetType  = 'Z';             % sheet with time x neurons matrix
distributionSelect = 'normal';% continuous outcome

% output folder for GLM results (per-file .xlsx)
pathOut = fullfile(dataPath, [scriptName]); % keeps results next to dataPath
if ~exist(pathOut, 'dir'), mkdir(pathOut); end

% ---------------- parameters that affect filenames ----------------
Nval = 10;
Tval=12000;
dtval = 0.05;
sigma_rval = 0.2;
tauval = 1.2;
wval = 0.25;
rho_targetval = 0.9;
smooth_sigmaval = 0.2;
seedval = 42;

% index for display
idx = 1;

% ---------------- iterate parameter grid ----------------
%for Tval = 2400:1200:12000
  for sparsityval = 0.05:0.05:0.3
    for snr_dbval = 0:5:10

      fileName = ['N_',num2str(Nval),'_T_',num2str(Tval),...
                  '_dt_',num2str(dtval),'_sparsity_',num2str(sparsityval),...
                  '_snr_',num2str(snr_dbval),'_sigma_',num2str(sigma_rval),...
                  '_tau_',num2str(tauval),'_w_',num2str(wval),'_rhoTarget_',...
                  num2str(rho_targetval),'_smoothSigma_',num2str(smooth_sigmaval),...
                  '_seed_',num2str(seedval),'.xlsx'];

      fprintf('Processing #%d -> %s\n', idx, fileName);
      infile = fullfile(pathIn, [dataType,'_', fileName]);
      if ~isfile(infile)
        warning('Missing file: %s. Skipping.', infile);
        idx = idx + 1;
        continue;
      end

      % read data: rows=time, cols=neurons expected on 'Z' sheet
      try
        data = readmatrix(infile, 'Sheet', sheetType);
      catch ME
        warning('Failed to read %s (sheet %s): %s. Skipping.', infile, sheetType, ME.message);
        idx = idx + 1;
        continue;
      end

      % convert to rows=time, cols=neurons if needed
      if size(data,1) < size(data,2)
        % your original did data' immediately. Keep original behaviour:
        data = data';
      end
      [numRows, numCols] = size(data);

      % ---- preallocate outputs (numCols x numCols) ----
      b1_CV1 = zeros(numCols, numCols);   b0_CV1 = zeros(1, numCols); lambda_CV1 = zeros(1, numCols);
      b1_CV2 = zeros(numCols, numCols);   b0_CV2 = zeros(1, numCols); lambda_CV2 = zeros(1, numCols);

      % new matrices for AIC / BIC / Custom
      b1_AIC = zeros(numCols, numCols);   b0_AIC = zeros(1, numCols); lambda_AIC = zeros(1, numCols);
      b1_BIC = zeros(numCols, numCols);   b0_BIC = zeros(1, numCols); lambda_BIC = zeros(1, numCols);
      b1_Custom = zeros(numCols, numCols); b0_Custom = zeros(1, numCols); lambda_Custom = zeros(1, numCols);

      mse_CV = nan(numCols, numCols);
      Rsqr_CV = nan(numCols, numCols);

      % create CV partition once (random K-fold)
      K = 10;
      cvp = cvpartition(numRows,'KFold',K);

      % ---- parfor across targets ----
      parfor neuron1 = 1:numCols
        y = data(:, neuron1);  % [T x 1] target

        % local buffers (to avoid slicing issues)
        local_bCV1 = zeros(numCols,1);
        local_b0CV1 = 0;
        local_lambdaCV1 = 0;

        local_bCV2 = zeros(numCols,1);
        local_b0CV2 = 0;
        local_lambdaCV2 = 0;

        local_bAIC = zeros(numCols,1);
        local_b0AIC = 0;
        local_lambdaAIC = 0;

        local_bBIC = zeros(numCols,1);
        local_b0BIC = 0;
        local_lambdaBIC = 0;

        local_bCustom = zeros(numCols,1);
        local_b0Custom = 0;
        local_lambdaCustom = 0;

        local_mse = nan(numCols,1);
        local_R2  = nan(numCols,1);

        for neuron2 = 1:numCols
          if neuron1 == neuron2
            % enforce zero self-weight
            local_bCV1(neuron2) = 0;
            local_bCV2(neuron2) = 0;
            local_bAIC(neuron2) = 0;
            local_bBIC(neuron2) = 0;
            local_bCustom(neuron2) = 0;
            continue;
          end

          X = data(:, neuron2);  % single-column predictor

          % Fit univariate LASSO GLM with CV
          try
            [bCV, statsCV] = lassoglm(X, y, distributionSelect, 'CV', cvp, 'Standardize', true);
          catch ME
            warning('lassoglm failed for target %d predictor %d: %s', neuron1, neuron2, ME.message);
            continue;
          end

          % sanity: statsCV fields existence
          if ~isfield(statsCV, 'IndexMinDeviance') || ~isfield(statsCV, 'Index1SE') || ~isfield(statsCV,'Lambda')
            warning('Unexpected statsCV structure for target %d predictor %d. Skipping pair.', neuron1, neuron2);
            continue;
          end

          % pick MinDeviance and 1SE as before
          idxMin = statsCV.IndexMinDeviance;
          idx1SE = statsCV.Index1SE;

          % coefficient matrix bCV is p x L (p==1)
          coef_min = bCV(:, idxMin);
          intercept_min = statsCV.Intercept(idxMin);
          lambda_min = statsCV.Lambda(idxMin);

          coef_1se = bCV(:, idx1SE);
          intercept_1se = statsCV.Intercept(idx1SE);
          lambda_1se = statsCV.Lambda(idx1SE);

          % store CV1 / CV2 results
          local_bCV1(neuron2) = coef_min;
          local_bCV2(neuron2) = coef_1se;

          % compute simple fit stats for this predictor-target at idxMin
          RSS = statsCV.Deviance(idxMin);
          TSS = sum((y - mean(y)).^2);
          local_R2(neuron2) = 1 - RSS/(TSS + eps);
          local_mse(neuron2) = RSS / max(1,numRows);

          local_b0CV1 = intercept_min;
          local_lambdaCV1 = lambda_min;
          local_b0CV2 = intercept_1se;
          local_lambdaCV2 = lambda_1se;

          % ---- compute AIC and BIC across lambdas for this pair ----
          dev = statsCV.Deviance; % numeric vector length L (may be RSS for gaussian)
          if isempty(dev) || all(isnan(dev))
            dev = inf(size(statsCV.Lambda));
          else
            dev(isnan(dev)) = inf;
          end

          if isfield(statsCV,'DF') && ~isempty(statsCV.DF)
            df = statsCV.DF;
            df(isnan(df)) = 0;
          else
            % fallback: estimate df per lambda by counting nonzero coeffs in bCV
            df = squeeze(sum(abs(bCV) > 1e-9, 1));
          end

          % AIC and BIC definitions used here: AIC = dev + 2*df ; BIC = dev + df*log(n)
          AICvals = dev + 2 .* df;
          BICvals = dev + df .* log(max(1,numRows));

          % pick lambda that minimizes AIC and BIC
          [~, idxAIC] = min(AICvals);
          [~, idxBIC] = min(BICvals);

          % save coefficients and intercepts at those indices
          local_bAIC(neuron2) = bCV(:, idxAIC);
          local_b0AIC = statsCV.Intercept(idxAIC);
          local_lambdaAIC = statsCV.Lambda(idxAIC);

          local_bBIC(neuron2) = bCV(:, idxBIC);
          local_b0BIC = statsCV.Intercept(idxBIC);
          local_lambdaBIC = statsCV.Lambda(idxBIC);

          % Custom rule: use BIC-minimizing lambda (change here if you want another rule)
          local_bCustom(neuron2) = local_bBIC(neuron2);
          local_b0Custom = local_b0BIC;
          local_lambdaCustom = local_lambdaBIC;

        end % neuron2

        % write local buffers into shared arrays (column = target)
        b1_CV1(:, neuron1) = local_bCV1;
        b0_CV1(:, neuron1) = local_b0CV1;
        lambda_CV1(:, neuron1) = local_lambdaCV1;

        b1_CV2(:, neuron1) = local_bCV2;
        b0_CV2(:, neuron1) = local_b0CV2;
        lambda_CV2(:, neuron1) = local_lambdaCV2;

        b1_AIC(:, neuron1) = local_bAIC;
        b0_AIC(:, neuron1) = local_b0AIC;
        lambda_AIC(:, neuron1) = local_lambdaAIC;

        b1_BIC(:, neuron1) = local_bBIC;
        b0_BIC(:, neuron1) = local_b0BIC;
        lambda_BIC(:, neuron1) = local_lambdaBIC;

        b1_Custom(:, neuron1) = local_bCustom;
        b0_Custom(:, neuron1) = local_b0Custom;
        lambda_Custom(:, neuron1) = local_lambdaCustom;

        mse_CV(:, neuron1) = local_mse;
        Rsqr_CV(:, neuron1) = local_R2;

      end % parfor neuron1

      % ensure diagonal zeros and no NaNs
      idxdiag = 1:(numCols+1):numCols^2;
      b1_CV1(idxdiag) = 0;
      b1_CV2(idxdiag) = 0;
      b1_AIC(idxdiag) = 0;
      b1_BIC(idxdiag) = 0;
      b1_Custom(idxdiag) = 0;

      b0_CV1(isnan(b0_CV1)) = 0;
      b0_CV2(isnan(b0_CV2)) = 0;
      b0_AIC(isnan(b0_AIC)) = 0;
      b0_BIC(isnan(b0_BIC)) = 0;
      b0_Custom(isnan(b0_Custom)) = 0;

      % ---------------- write outputs ----------------
      outfile = fullfile(pathOut, [scriptName,'_', fileName]);

      % CV1/CV2
      writematrix(b1_CV1, outfile, 'Sheet', 'b1_CV1');
      writematrix(b0_CV1, outfile, 'Sheet', 'b0_CV1');
      writematrix(lambda_CV1, outfile, 'Sheet', 'lambda_CV1');

      writematrix(b1_CV2, outfile, 'Sheet', 'b1_CV2');
      writematrix(b0_CV2, outfile, 'Sheet', 'b0_CV2');
      writematrix(lambda_CV2, outfile, 'Sheet', 'lambda_CV2');

      % AIC
      writematrix(b1_AIC, outfile, 'Sheet', 'b1_AIC');
      writematrix(b0_AIC, outfile, 'Sheet', 'b0_AIC');
      writematrix(lambda_AIC, outfile, 'Sheet', 'lambda_AIC');

      % BIC
      writematrix(b1_BIC, outfile, 'Sheet', 'b1_BIC');
      writematrix(b0_BIC, outfile, 'Sheet', 'b0_BIC');
      writematrix(lambda_BIC, outfile, 'Sheet', 'lambda_BIC');

      % Custom
      writematrix(b1_Custom, outfile, 'Sheet', 'b1_Custom');
      writematrix(b0_Custom, outfile, 'Sheet', 'b0_Custom');
      writematrix(lambda_Custom, outfile, 'Sheet', 'lambda_Custom');

      % extras
      writematrix(mse_CV, outfile, 'Sheet', 'mse_CV');
      writematrix(Rsqr_CV, outfile, 'Sheet', 'Rsqr_CV');

      fprintf('Saved: %s\n', outfile);
      idx = idx + 1;
    end % snr_dbval
  end % sparsity
%end % Tval

fprintf('All done. Results written to %s\n', pathOut);