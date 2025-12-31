%% fitglm_allpredictors_network.m
clear; close all; clc; warning('off');
scriptName = mfilename;

%% Ensure workers can see current code path (parfor safety)
p = gcp('nocreate');
if isempty(p); p = parpool; end
% Replicate current path on workers
pctRunOnAll addpath(genpath(pwd));

%% ---- Inputs (same dialog as before)
pathIn = ['..',filesep,'data',filesep,'groundTruth',filesep];

pathOut = ['..',filesep,'data',filesep,'spike',filesep,scriptName,filesep];
if ~exist(pathOut , 'dir')
    mkdir(pathOut );
end
spikeSheet='rate_smooth';
distributionSelect='poisson';% spikes must have poission or binomial

Nval=10;%10-30
Tval=2400;%2400–12 000
dtval=0.05;
sparsityval=0.05;%0.05 – 0.30
snr_dbval=0;%0–10 dB ;0 dB = very noisy, 5 dB = medium, 10 dB = high. 
seedval=42;
wval=0.25;
sigma_rval=0.2;
rho_targetval=0.9;
tauval=1.2;
smooth_sigmaval=0.2;

%for Nval=10:10:30;
    Nval=30;
    for Tval=2400:1200:12000
        for sparsityval=0.05:0.05:0.3
            for snr_dbval=0:5:10

 fileName=['groundTruth_N_',num2str(Nval),'_T_',num2str(Tval),...
                        '_dt_',num2str(dtval),'_sparsity_',num2str(sparsityval),...
                        '_snr_',num2str(snr_dbval),'_sigma_',num2str(sigma_rval),...
                        '_tau_',num2str(tauval),'_w_',num2str(wval),'_rhoTarget_',...
                        num2str(rho_targetval),'_smoothSigma_',num2str(smooth_sigmaval),...
                        '_seed_',num2str(seedval),'.xlsx'];

data = readmatrix([pathIn,fileName,],'Sheet',spikeSheet);
data=data';%columns are z score traces for individual neurons
dt=dtval;

dataCa = data(:,1:end);%Ca fluoresence z-scores

[numRows, numCols] = size(dataCa);

%% ---- Preallocate all outputs (n×n matrices, same names/sheets as your original)
b1_fitglm       = NaN(numCols, numCols);   % slope (beta1: source->target)
b0_fitglm       = NaN(numCols, numCols);   % intercept (beta0)
se1_fitglm      = NaN(numCols, numCols);   % SE(beta1)
se0_fitglm      = NaN(numCols, numCols);   % SE(beta0)
t1_fitglm       = NaN(numCols, numCols);   % t-stat beta1
t0_fitglm       = NaN(numCols, numCols);   % t-stat beta0
p1_fitglm       = NaN(numCols, numCols);   % p-value beta1
p0_fitglm       = NaN(numCols, numCols);   % p-value beta0
deviance_fitglm = NaN(numCols, numCols);   % SSE for 'normal'
df_fitglm       = NaN(numCols, numCols);   % residual DOF
AIC_fitglm      = NaN(numCols, numCols);
BIC_fitglm      = NaN(numCols, numCols);

% Cross-validated / holdout metrics (replicated by column to keep shape)
mse_CV10  = NaN(numCols, numCols);
R2_CV10   = NaN(numCols, numCols);
mse_HO20  = NaN(numCols, numCols);
R2_HO20   = NaN(numCols, numCols);

%% ---- Settings
K = 10;                    % K-fold for CV
STANDARDIZE_FOR_CV = true; % standardize predictors only in CV/holdout (train stats)

%% ---- Multi-predictor fits (target = neuron1, predictors = all others)
% parpool('threads'); % optional

parfor neuron1 = 1:numCols
    % ---------------- Response and predictors ----------------
    y = dataCa(:, neuron1);
    predIdx = setdiff(1:numCols, neuron1);   % use all other columns as predictors
    X = dataCa(:, predIdx);

    % ------------------ Fit GLM (unregularized, on original scale) ------------------
    % Intercept is included by default; 'Intercept',true is explicit and valid.
    mdl = fitglm(X, y, 'Distribution','normal', 'Intercept', true);

    % Coefficients table: Row 1 = Intercept, rows 2:end follow column order of X (predIdx)
    coefs = mdl.Coefficients;

    % Map estimates and stats back to global n×n matrices
    % Initialize local column vectors for this target
    b1_col  = NaN(numCols,1);  b0_col  = NaN(numCols,1);
    se1_col = NaN(numCols,1);  se0_col = NaN(numCols,1);
    t1_col  = NaN(numCols,1);  t0_col  = NaN(numCols,1);
    p1_col  = NaN(numCols,1);  p0_col  = NaN(numCols,1);

    % Intercept: replicate down the column, but NaN on the diagonal (no self-coupling)
    b0_scalar  = coefs.Estimate(1);
    se0_scalar = coefs.SE(1);
    t0_scalar  = coefs.tStat(1);
    p0_scalar  = coefs.pValue(1);

    b0_col(:)  = b0_scalar;  b0_col(neuron1)  = NaN;
    se0_col(:) = se0_scalar; se0_col(neuron1) = NaN;
    t0_col(:)  = t0_scalar;  t0_col(neuron1)  = NaN;
    p0_col(:)  = p0_scalar;  p0_col(neuron1)  = NaN;

    % Slopes and their SE/t/p for the (n-1) predictors
    % Rows 2..end in 'coefs' correspond to X(:, predIdx) in that order
    b1_vals  = coefs.Estimate(2:end);
    se1_vals = coefs.SE(2:end);
    t1_vals  = coefs.tStat(2:end);
    p1_vals  = coefs.pValue(2:end);

    % Place them into the correct predictor rows
    b1_col(predIdx)  = b1_vals;
    se1_col(predIdx) = se1_vals;
    t1_col(predIdx)  = t1_vals;
    p1_col(predIdx)  = p1_vals;

    % ------------------ Model-level metrics (replicated) ----------------------
    dev_scalar = mdl.Deviance;        % SSE for 'normal'
    df_scalar  = mdl.DFE;
    aic_scalar = mdl.ModelCriterion.AIC;
    bic_scalar = mdl.ModelCriterion.BIC;

    dev_col = repmat(dev_scalar, numCols, 1);
    df_col  = repmat(df_scalar,  numCols, 1);
    aic_col = repmat(aic_scalar, numCols, 1);
    bic_col = repmat(bic_scalar, numCols, 1);

    % Optional: set diagonal NaN for cosmetic consistency
    dev_col(neuron1) = NaN;
    df_col(neuron1)  = NaN;
    aic_col(neuron1) = NaN;
    bic_col(neuron1) = NaN;

    % ------------------ 10-fold CV (multi-predictor) -------------------
    c = cvpartition(numRows, 'KFold', K);
    yhat_all = NaN(numRows,1);
    for k = 1:K
        tr = training(c,k); te = test(c,k);
        if STANDARDIZE_FOR_CV
            % Train-only mean/std
            mu  = mean(X(tr,:), 1);
            sig = std(X(tr,:), 0, 1);
            sig(sig==0) = 1; % avoid divide-by-zero
            XtrZ = (X(tr,:) - mu) ./ sig;
            XteZ = (X(te,:) - mu) ./ sig;
            mdl_k   = fitglm(XtrZ, y(tr), 'Distribution','normal', 'Intercept', true);
            yhat_all(te) = predict(mdl_k, XteZ);
        else
            mdl_k   = fitglm(X(tr,:), y(tr), 'Distribution','normal', 'Intercept', true);
            yhat_all(te) = predict(mdl_k, X(te,:));
        end
    end
    mseCV_scalar = mean((yhat_all - y).^2, 'omitnan');
    ymean = mean(y);
    TSS   = sum((y - ymean).^2);
    R2CV_scalar = 1 - (mseCV_scalar * numRows) / TSS;

    mseCV_col = repmat(mseCV_scalar, numCols, 1);
    R2CV_col  = repmat(R2CV_scalar,  numCols, 1);
    mseCV_col(neuron1) = NaN;
    R2CV_col(neuron1)  = NaN;

    % ------------------ Hold-out 20% (multi-predictor) -------------------
    ho = cvpartition(numRows, 'HoldOut', 0.2);
    tr = training(ho); te = test(ho);
    if STANDARDIZE_FOR_CV
        muH  = mean(X(tr,:), 1);
        sigH = std(X(tr,:), 0, 1);
        sigH(sigH==0) = 1;
        XtrZ = (X(tr,:) - muH) ./ sigH;
        XteZ = (X(te,:) - muH) ./ sigH;
        mdl_ho  = fitglm(XtrZ, y(tr), 'Distribution','normal', 'Intercept', true);
        yhat_te = predict(mdl_ho, XteZ);
    else
        mdl_ho  = fitglm(X(tr,:), y(tr), 'Distribution','normal', 'Intercept', true);
        yhat_te = predict(mdl_ho, X(te,:));
    end
    mseHO_scalar = mean((yhat_te - y(te)).^2);
    ymean_te = mean(y(te));
    TSS_te   = sum((y(te) - ymean_te).^2);
    R2HO_scalar = 1 - (mseHO_scalar * sum(te)) / TSS_te;

    mseHO_col = repmat(mseHO_scalar, numCols, 1);
    R2HO_col  = repmat(R2HO_scalar,  numCols, 1);
    mseHO_col(neuron1) = NaN;
    R2HO_col(neuron1)  = NaN;

    % ----------------- Commit this target column -------------------------
    b1_fitglm(:,neuron1)       = b1_col;
    b0_fitglm(:,neuron1)       = b0_col;
    se1_fitglm(:,neuron1)      = se1_col;
    se0_fitglm(:,neuron1)      = se0_col;
    t1_fitglm(:,neuron1)       = t1_col;
    t0_fitglm(:,neuron1)       = t0_col;
    p1_fitglm(:,neuron1)       = p1_col;
    p0_fitglm(:,neuron1)       = p0_col;

    deviance_fitglm(:,neuron1) = dev_col;
    df_fitglm(:,neuron1)       = df_col;
    AIC_fitglm(:,neuron1)      = aic_col;
    BIC_fitglm(:,neuron1)      = bic_col;

    mse_CV10(:,neuron1)        = mseCV_col;
    R2_CV10(:,neuron1)         = R2CV_col;
    mse_HO20(:,neuron1)        = mseHO_col;
    R2_HO20(:,neuron1)         = R2HO_col;
end

%% ---- Save to Excel (same sheet names)
fileOut=[pathOut,extractAfter(fileName,'groundTruth')];
writematrix(b1_fitglm,       fileOut,'Sheet','beta1');
writematrix(b0_fitglm,       fileOut,'Sheet','beta0');
writematrix(se1_fitglm,      fileOut,'Sheet','se_beta1');
writematrix(se0_fitglm,      fileOut,'Sheet','se_beta0');
writematrix(t1_fitglm,       fileOut,'Sheet','t_beta1');
writematrix(t0_fitglm,       fileOut,'Sheet','t_beta0');
writematrix(p1_fitglm,       fileOut,'Sheet','p_beta1');
writematrix(p0_fitglm,       fileOut,'Sheet','p_beta0');
writematrix(deviance_fitglm, fileOut,'Sheet','Deviance');
writematrix(df_fitglm,       fileOut,'Sheet','DFE');
writematrix(AIC_fitglm,      fileOut,'Sheet','AIC');
writematrix(BIC_fitglm,      fileOut,'Sheet','BIC');
writematrix(mse_CV10,        fileOut,'Sheet','MSE_CV10');
writematrix(R2_CV10,         fileOut,'Sheet','R2_CV10');
writematrix(mse_HO20,        fileOut,'Sheet','MSE_Holdout');
writematrix(R2_HO20,         fileOut,'Sheet','R2_Holdout');
            end
        end
    end
    %end