%take ground truth synthetic data and run them through GLM to get the
%estimated connectivity
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

%% Selective Data Formatting
pathIn = ['..',filesep,'data',filesep,'groundTruth',filesep];

pathOut = ['..',filesep,'data',filesep,scriptName,filesep];
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

%for Nval=10:10:30
    Nval=30;
    for Tval=2400:1200:12000
        %Tval=2400;
        for sparsityval=0.05:0.05:0.3
            %sparsityval=0.05;
            for snr_dbval=0:5:10
                %snr_dbval=0;
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

                %% General Linear Model
                % See EIER for basic linear regression model

                    % clear bestCV_1 interceptCV_1  bestLambdaDev bestCV_2 interceptCV_2
                    % clear bestLambda1SE bestMSE R2_CV bestB_AIC intercept_AIC bestB_BIC
                    % clear intercept_BIC bestB_custom intercept_custom mseCustomMin R2_Custom;

                parfor neuron1=1:numCols
                    y = dataCa(:,neuron1); %select the neuron to predict its activity
                    % Local buffers for parfor slicing
                    bestCV_1= NaN(numCols, 1);
                    interceptCV_1= NaN(numCols, 1);
                    bestLambdaDev= NaN(numCols, 1);
                    bestCV_2= NaN(numCols, 1);
                    interceptCV_2= NaN(numCols, 1);
                    bestLambda1SE= NaN(numCols, 1);
                    bestMSE= NaN(numCols, 1);
                    R2_CV= NaN(numCols, 1);
                    bestB_AIC= NaN(numCols, 1);
                    intercept_AIC= NaN(numCols, 1);
                    bestB_BIC= NaN(numCols, 1);
                    intercept_BIC= NaN(numCols, 1);
                    bestB_custom= NaN(numCols, 1);
                    intercept_custom= NaN(numCols, 1);
                    mseCustomMin= NaN(numCols, 1);
                    R2_Custom= NaN(numCols, 1);

                    for neuron2=1:numCols
                        if neuron1 == neuron2, continue; end
                        X = dataCa(:,neuron2);

                        % 10-fold cross validation CV; change 'normal' to 'binomial' etc. if needed
                        [bCV, statsCV] = lassoglm(X, y,distributionSelect, 'CV', 10);  
                        bestLambdaDev(neuron2)= statsCV.LambdaMinDeviance;
                        bestLambda1SE(neuron2)=statsCV.Lambda1SE;
                        %Methiod 1: Get the best Lambda (minimum deviance) 
                        idxLambdaMinDevianceCV =statsCV.IndexMinDeviance;
                        bestCV_1(neuron2) = bCV(:, idxLambdaMinDevianceCV);  % The best coefficients
                        %The intercept is stored separately:
                        interceptCV_1(neuron2) = statsCV.Intercept(idxLambdaMinDevianceCV);

                        %Method 2: Regularized model within one standard error to avoid overfitting:
                        idxLambda1SECV = statsCV.Index1SE;
                        bestCV_2(neuron2) = bCV(:, idxLambda1SECV);
                        interceptCV_2(neuron2) = statsCV.Intercept(idxLambda1SECV);
                        %So your full prediction model is:
                        %yhat = X * bestB + intercept;
                        %For classification (e.g. binomial), you'd apply a link function:
                        %yhatProb = 1 ./ (1 + exp(-(X * bestB + intercept)));  % sigmoid

                        % Statistics on fiting
                        mseCV = statsCV.Deviance / numRows;  % Deviance is SSE here
                        bestMSE(neuron2) = mseCV(statsCV.IndexMinDeviance);
                        %You can calculate pseudo-R² for regression:
                        ymean = mean(y);
                        TSS = sum((y - ymean).^2);
                        RSS = statsCV.Deviance(idxLambdaMinDevianceCV);  % from statsCV.Deviance(bestIdx)
                        R2_CV(neuron2) = 1 - RSS / TSS;

                        %%%%%%% AIC and BIC (Akaike & Bayesian Information Criteria)
                        %%%%%%%% Both AIC and BIC balance goodness of fit and model complexity.
                        AIC = statsCV.Deviance + 2 * statsCV.DF;
                        BIC = statsCV.Deviance + log(numRows) * statsCV.DF;

                        [minAIC, idxAIC] = min(AIC);
                        [minBIC, idxBIC] = min(BIC);

                        bestB_AIC(neuron2) = bCV(:, idxAIC);
                        intercept_AIC(neuron2) = statsCV.Intercept(idxAIC);

                        bestB_BIC(neuron2) = bCV(:, idxBIC);
                        intercept_BIC(neuron2) = statsCV.Intercept(idxBIC);   

                        %%%%%%%%%%%Custom Validation (e.g. holdout or K-fold)
                        %%%%%%%%%%%%Example with Holdout Validation:
                        % Split data
                        customValid = cvpartition(size(X,1), 'HoldOut', 0.2);
                        Xtrain = X(training(customValid), :);
                        ytrain = y(training(customValid));
                        Xtest = X(test(customValid), :);
                        ytest = y(test(customValid));
                        [bCustom, statsCustom] = lassoglm(Xtrain, ytrain, distributionSelect, 'Lambda', logspace(-4, 1, 100));% Fit model
                        % Split the data
                        % Predict on test set
                        preds = Xtest * bCustom + statsCustom.Intercept;  % size: [nTest x numLambda]

                        % Compute mean squared error for each lambda
                        mseCustom = mean((preds - ytest).^2);

                        % Choose lambda with lowest MSE
                        [mseCustomMin(neuron2), IdxCustom] = min(mseCustom);
                        bestB_custom(neuron2) = bCustom(:, IdxCustom);
                        intercept_custom(neuron2) = statsCustom.Intercept(IdxCustom);

                        %You can calculate pseudo-R² for regression:
                         ymean_test = mean(ytest);
                        TSS = sum((ytest - ymean_test).^2);
                        RSS = mseCustom(IdxCustom) * numel(ytest);
                        R2_Custom(neuron2) = 1 - RSS / TSS;
                    end;%neuron 2 loop
                    %         % Assign back into global result matrices
                   b1_CV1(:, neuron1)= bestCV_1;
                   b0_CV1(:, neuron1)= interceptCV_1;
                   lambda_CV1(:, neuron1)=bestLambdaDev;
                   b1_CV2(:, neuron1)=bestCV_2;
                   b0_CV2(:, neuron1)= interceptCV_2;
                   lambda_CV2(:, neuron1)= bestLambda1SE;
                   mse_CV(:, neuron1)=bestMSE;
                   Rsqr_CV(:, neuron1)= R2_CV;
                   b1_AIC(:, neuron1)= bestB_AIC;
                   b0_AIC(:, neuron1)= intercept_AIC;
                   b1_BIC(:, neuron1)=bestB_BIC;
                   b0_BIC(:, neuron1)= intercept_BIC;
                   b1_Custom(:, neuron1)= bestB_custom;
                   b0_Custom(:, neuron1)=intercept_custom;
                   mse_Custom(:, neuron1)= mseCustomMin;
                   Rsqr_Custom(:, neuron1)=R2_Custom;
                end; % end all neuron1

                %Each column is the beta coefficients for neuron1
                fileOut=[pathOut,scriptName,extractAfter(fileName,'groundTruth')];
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
                writematrix( b1_Custom,fileOut,'Sheet','b1_Custom')
                writematrix(b0_Custom,fileOut,'Sheet','b0_Custom')
                writematrix(mse_Custom,fileOut,'Sheet','mse_Custom')
                writematrix(Rsqr_Custom,fileOut,'Sheet','Rsqr_Custom')

            end
        end
    end
%end; %Nval