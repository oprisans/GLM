clear all;
close all;
clc;
warning('off')
scriptName = mfilename;

%% Ensure workers can see current code path (parfor safety)
p = gcp('nocreate');
if isempty(p); p = parpool; end
% Replicate current path on workers
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
pathOut = fullfile(dataPath,scriptName); % keeps results next to dataPath
if ~exist(pathOut, 'dir'), mkdir(pathOut); end

Nval=10;%10-30
Tval=12000;%2400–12 000
dtval=0.05;
sparsityval=0.05;%0.05 – 0.30
snr_dbval=0;%0–10 dB ;0 dB = very noisy, 5 dB = medium, 10 dB = high. 
seedval=42;
wval=0.25;
sigma_rval=0.2;
rho_targetval=0.9;
tauval=1.2;
smooth_sigmaval=0.2;

for Nval=10:10:30;
   %for Tval=2400:1200:12000
        for sparsityval=0.05:0.05:0.5
            for snr_dbval=0:5:10

 fileName=['groundTruth_N_',num2str(Nval),'_T_',num2str(Tval),...
                        '_dt_',num2str(dtval),'_sparsity_',num2str(sparsityval),...
                        '_snr_',num2str(snr_dbval),'_sigma_',num2str(sigma_rval),...
                        '_tau_',num2str(tauval),'_w_',num2str(wval),'_rhoTarget_',...
                        num2str(rho_targetval),'_smoothSigma_',num2str(smooth_sigmaval),...
                        '_seed_',num2str(seedval),'.xlsx'];

data = readmatrix(fullfile(pathIn,fileName),'Sheet',sheetType);
data=data';%columns are z score traces for individual neurons
dt=dtval;

dataCa = data(:,1:end);%Ca fluoresence z-scores
[numRows, numCols] = size(dataCa);

%% General Linear Model (univariate response; all other columns as predictors)
% Continuous data -> 'normal'

% --------- Preallocate global result matrices (n x n) ----------
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

parfor neuron1 = 1:numCols
    % ---------------- Response and predictors ----------------
    y = dataCa(:, neuron1);
    predIdx = setdiff(1:numCols, neuron1); % use all other columns as predictors
    X = dataCa(:, predIdx);

    % Local buffers mapped to predictor rows; will be assigned to column neuron1
    bestCV_1             = NaN(numCols,1);
    interceptCV_1        = NaN(numCols,1);
    bestLambdaDev_col    = NaN(numCols,1);

    bestCV_2             = NaN(numCols,1);
    interceptCV_2        = NaN(numCols,1);
    bestLambda1SE_col    = NaN(numCols,1);

    bestMSE_col          = NaN(numCols,1);
    R2_CV_col            = NaN(numCols,1);

    bestB_AIC_col        = NaN(numCols,1);
    intercept_AIC_col    = NaN(numCols,1);

    bestB_BIC_col        = NaN(numCols,1);
    intercept_BIC_col    = NaN(numCols,1);

    bestB_custom_col     = NaN(numCols,1);
    intercept_custom_col = NaN(numCols,1);
    mseCustomMin_col     = NaN(numCols,1);
    R2_Custom_col        = NaN(numCols,1);

    % ===================== 10-fold Cross-Validation ======================
    [bCV, statsCV] = lassoglm(X, y, distributionSelect, 'CV', 10,'Standardize',true);
    idxMinDev = statsCV.IndexMinDeviance;
    idx1SE    = statsCV.Index1SE;

    % Map coefficients (n-1) -> length-n with NaN at self index
    coefMin = NaN(numCols,1);  coefMin(predIdx) = bCV(:, idxMinDev);
    coef1SE = NaN(numCols,1);  coef1SE(predIdx) = bCV(:, idx1SE);

    bestCV_1 = coefMin;
    bestCV_2 = coef1SE;

    % Intercepts: replicate down column; set diagonal to NaN (no self-coupling)
    interceptCV_1 = repmat(statsCV.Intercept(idxMinDev), numCols, 1);
    interceptCV_1(neuron1) = NaN;

    interceptCV_2 = repmat(statsCV.Intercept(idx1SE), numCols, 1);
    interceptCV_2(neuron1) = NaN;

    % Lambdas: replicate; diagonal NaN for cosmetic consistency
    bestLambdaDev_col = repmat(statsCV.LambdaMinDeviance, numCols, 1);
    bestLambdaDev_col(neuron1) = NaN;

    bestLambda1SE_col = repmat(statsCV.Lambda1SE, numCols, 1);
    bestLambda1SE_col(neuron1) = NaN;

    % CV Deviance -> MSE; then R² vs in-sample TSS (matches your original)
    mse_all = statsCV.Deviance ./ numRows;
    bestMSE_val = mse_all(idxMinDev);
    bestMSE_col = repmat(bestMSE_val, numCols, 1);
    bestMSE_col(neuron1) = NaN;

    ymean = mean(y);
    TSS = sum((y - ymean).^2);
    RSS = statsCV.Deviance(idxMinDev);
    R2_val = 1 - RSS / TSS;
    R2_CV_col = repmat(R2_val, numCols, 1);
    R2_CV_col(neuron1) = NaN;

    % =========================== AIC and BIC =============================
    % Compute AIC/BIC from CV path quantities; pick best lambda
    AIC = statsCV.Deviance + 2 * statsCV.DF;
    BIC = statsCV.Deviance + log(numRows) * statsCV.DF;

    [~, idxAIC] = min(AIC);
    [~, idxBIC] = min(BIC);

    coefAIC = NaN(numCols,1);   coefAIC(predIdx) = bCV(:, idxAIC);
    coefBIC = NaN(numCols,1);   coefBIC(predIdx) = bCV(:, idxBIC);

    bestB_AIC_col = coefAIC;
    bestB_BIC_col = coefBIC;

    intercept_AIC_col = repmat(statsCV.Intercept(idxAIC), numCols, 1);
    intercept_AIC_col(neuron1) = NaN;

    intercept_BIC_col = repmat(statsCV.Intercept(idxBIC), numCols, 1);
    intercept_BIC_col(neuron1) = NaN;

    % ========================= Custom Holdout ============================
    customValid = cvpartition(size(X,1), 'HoldOut', 0.2);
    Xtrain = X(training(customValid), :);
    ytrain = y(training(customValid));
    Xtest  = X(test(customValid), :);
    ytest  = y(test(customValid));

    [bCustom, statsCustom] = lassoglm(Xtrain, ytrain, distributionSelect, ...
                                      'Lambda', logspace(-4, 1, 100),'Standardize',true);

    % Predictions for all lambdas (use bsxfun for wide compatibility)
    % preds size: [nTest x nLambda]
interceptRow = reshape(statsCustom.Intercept, 1, []);
preds        = Xtest * bCustom + interceptRow;   % implicit broadcasting    

    mseCustom = mean((preds - ytest).^2, 1);
    [mseMin, idxCustom] = min(mseCustom);

    coefCustom = NaN(numCols,1);
    coefCustom(predIdx) = bCustom(:, idxCustom);

    bestB_custom_col = coefCustom;

    intercept_custom_col = repmat(statsCustom.Intercept(idxCustom), numCols, 1);
    intercept_custom_col(neuron1) = NaN;

    mseCustomMin_col = repmat(mseMin, numCols, 1);
    mseCustomMin_col(neuron1) = NaN;

    ymean_test = mean(ytest);
    TSS_test = sum((ytest - ymean_test).^2);
    RSS_test = mseMin * numel(ytest);
    R2_Custom_val = 1 - RSS_test / TSS_test;
    R2_Custom_col = repmat(R2_Custom_val, numCols, 1);
    R2_Custom_col(neuron1) = NaN;

    % ----------------- Assign back into global result matrices -----------
    b1_CV1(:, neuron1)      = bestCV_1;
    b0_CV1(:, neuron1)      = interceptCV_1;
    lambda_CV1(:, neuron1)  = bestLambdaDev_col;

    b1_CV2(:, neuron1)      = bestCV_2;
    b0_CV2(:, neuron1)      = interceptCV_2;
    lambda_CV2(:, neuron1)  = bestLambda1SE_col;

    mse_CV(:, neuron1)      = bestMSE_col;
    Rsqr_CV(:, neuron1)     = R2_CV_col;

    b1_AIC(:, neuron1)      = bestB_AIC_col;
    b0_AIC(:, neuron1)      = intercept_AIC_col;

    b1_BIC(:, neuron1)      = bestB_BIC_col;
    b0_BIC(:, neuron1)      = intercept_BIC_col;

    b1_Custom(:, neuron1)   = bestB_custom_col;
    b0_Custom(:, neuron1)   = intercept_custom_col;
    mse_Custom(:, neuron1)  = mseCustomMin_col;
    Rsqr_Custom(:, neuron1) = R2_Custom_col;
end % parfor neuron1

% Each column is the set of coefficients for 'neuron1' as response
fileOut=fullfile(pathOut, [scriptName,'_', fileName]);
writematrix(b1_CV1,fileOut,'Sheet','b1_CV1')
writematrix(b0_CV1,fileOut,'Sheet','b0_CV')
writematrix(lambda_CV1,fileOut,'Sheet','lambda_CV1')
writematrix(b1_CV2,fileOut,'Sheet','b1_CV2')
writematrix(b0_CV2,fileOut,'Sheet','b0_CV2')
writematrix(lambda_CV2,fileOut,'Sheet','lambda_CV2')
writematrix(mse_CV,fileOut,'Sheet','mse_CV')
writematrix(Rsqr_CV,fileOut,'Sheet','Rsqr_CV')
writematrix(b1_AIC,fileOut,'Sheet','b1_AIC')
writematrix(b0_AIC,fileOut,'Sheet','b0_AIC')
writematrix(b1_BIC,fileOut,'Sheet','b1_BIC')
writematrix(b0_BIC,fileOut,'Sheet','b0_BIC')
writematrix(b1_Custom,fileOut,'Sheet','b1_Custom')
writematrix(b0_Custom,fileOut,'Sheet','b0_Custom')
writematrix(mse_Custom,fileOut,'Sheet','mse_Custom')
writematrix(Rsqr_Custom,fileOut,'Sheet','Rsqr_Custom')

            end
        end
    end
    %end
