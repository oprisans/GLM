%% fitglm_pairwise_network.m
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

%% ---- Preallocate all outputs (n×n matrices, like your original)
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

% Cross-validated metrics
mse_CV10  = NaN(numCols, numCols);
R2_CV10   = NaN(numCols, numCols);
mse_HO20  = NaN(numCols, numCols);
R2_HO20   = NaN(numCols, numCols);

%% ---- Parallel pairwise fits (target = neuron1, predictor = neuron2)
% Tip: make sure you have a parallel pool if you want speedup
% parpool('threads');  % optional
K = 10; % k-fold for CV
parfor neuron1 = 1:numCols
    y = dataCa(:, neuron1);

    % local buffers per target (parfor-friendly)
    b1_loc  = NaN(numCols,1);  b0_loc  = NaN(numCols,1);
    se1_loc = NaN(numCols,1);  se0_loc = NaN(numCols,1);
    t1_loc  = NaN(numCols,1);  t0_loc  = NaN(numCols,1);
    p1_loc  = NaN(numCols,1);  p0_loc  = NaN(numCols,1);
    dev_loc = NaN(numCols,1);  df_loc  = NaN(numCols,1);
    aic_loc = NaN(numCols,1);  bic_loc = NaN(numCols,1);
    mseCV_loc = NaN(numCols,1); R2CV_loc = NaN(numCols,1);
    mseHO_loc = NaN(numCols,1); R2HO_loc = NaN(numCols,1);

    for neuron2 = 1:numCols
        if neuron1 == neuron2, continue; end

        X = dataCa(:, neuron2);  % single-predictor univariate GLM

        % ------------------ Fit GLM (unregularized) ------------------
        mdl = fitglm(X, y, 'Distribution',distributionSelect, 'Intercept', true);

        % coefficients (rows: Intercept; x1)
        coefs = mdl.Coefficients;
        b0_loc(neuron2)  = coefs.Estimate(1);
        b1_loc(neuron2)  = coefs.Estimate(2);
        se0_loc(neuron2) = coefs.SE(1);
        se1_loc(neuron2) = coefs.SE(2);
        t0_loc(neuron2)  = coefs.tStat(1);
        t1_loc(neuron2)  = coefs.tStat(2);
        p0_loc(neuron2)  = coefs.pValue(1);
        p1_loc(neuron2)  = coefs.pValue(2);

        dev_loc(neuron2) = mdl.Deviance;            % SSE for normal
        df_loc(neuron2)  = mdl.DFE;
        aic_loc(neuron2) = mdl.ModelCriterion.AIC;
        bic_loc(neuron2) = mdl.ModelCriterion.BIC;

        % ------------------ 10-fold CV (manual) ----------------------
        c = cvpartition(numRows, 'KFold', K);
        yhat_all = NaN(numRows,1);
        for k = 1:K
            tr = training(c,k); te = test(c,k);
            mdl_k = fitglm(X(tr), y(tr), 'Distribution',distributionSelect, 'Intercept', true);
            yhat_all(te) = predict(mdl_k, X(te));
        end
        mseCV = mean((yhat_all - y).^2, 'omitnan');
        mseCV_loc(neuron2) = mseCV;

        ymean = mean(y);
        TSS   = sum((y - ymean).^2);
        R2CV_loc(neuron2) = 1 - (mseCV * numRows) / TSS;

        % ------------------ Hold-out 20% -----------------------------
        ho = cvpartition(numRows, 'HoldOut', 0.2);
        tr = training(ho); te = test(ho);
        mdl_ho = fitglm(X(tr), y(tr), 'Distribution',distributionSelect, 'Intercept', true);
        yhat_te = predict(mdl_ho, X(te));
        mseHO = mean((yhat_te - y(te)).^2);
        mseHO_loc(neuron2) = mseHO;

        ymean_te = mean(y(te));
        TSS_te   = sum((y(te) - ymean_te).^2);
        R2HO_loc(neuron2) = 1 - (mseHO * sum(te)) / TSS_te;
    end

    % commit this target column
    b1_fitglm(:,neuron1)       = b1_loc;
    b0_fitglm(:,neuron1)       = b0_loc;
    se1_fitglm(:,neuron1)      = se1_loc;
    se0_fitglm(:,neuron1)      = se0_loc;
    t1_fitglm(:,neuron1)       = t1_loc;
    t0_fitglm(:,neuron1)       = t0_loc;
    p1_fitglm(:,neuron1)       = p1_loc;
    p0_fitglm(:,neuron1)       = p0_loc;
    deviance_fitglm(:,neuron1) = dev_loc;
    df_fitglm(:,neuron1)       = df_loc;
    AIC_fitglm(:,neuron1)      = aic_loc;
    BIC_fitglm(:,neuron1)      = bic_loc;
    mse_CV10(:,neuron1)        = mseCV_loc;
    R2_CV10(:,neuron1)         = R2CV_loc;
    mse_HO20(:,neuron1)        = mseHO_loc;
    R2_HO20(:,neuron1)         = R2HO_loc;
end

%% ---- Save to Excel
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

%disp('Done: results written to:'); disp(fileOut);

            end
        end
    end
    %end