%This simulator matches the benchmarking aim: you know the ground-truth directed matrix A, 
% you pass Z and rate_smooth through your pairwise/multi-predictor 
% linear models (with/without L1), and you score recovery (AUPRC/ROC, 
% weight correlation, SHD).

% Get the current script name
scriptName = mfilename('fullpath');
% Extract the directory path
scriptPath = fileparts(scriptName);
% Move one step up in the directory hierarchy
parentPath = fileparts(scriptPath);
% Create the new path to the 'data' directory
dataPath = fullfile(parentPath, 'data');
resultsPath=fullfile(parentPath, 'results');

dataType='groundTruth';
pathOut = fullfile(dataPath,dataType);
if ~exist(pathOut , 'dir')
    mkdir(pathOut );
end

Nval=10;%10-30
Tval=2400;%2400–12 000
dtval=0.05;
sparsityval=0.05;%0.05 – 0.30
snr_sbval=0;%0–10 dB ;0 dB = very noisy, 5 dB = medium, 10 dB = high. 
seedval=42;
wval=0.25;
sigma_rval=0.2;
rho_targetval=0.9;
tauval=1.2;
smooth_sigmaval=0.2;

Nval=10;
for Tval = 2400:1200:12000
  for sparsityval = 0.05:0.05:0.3
    for snr_dbval = 0:5:10

                    sim = simulate_var_calcium('N',Nval,'T',Tval,'dt',dtval,'sparsity',sparsityval,...
                        'snr_db',snr_dbval,'sigma_r',sigma_rval,'tau',tauval,'w',wval,'rho_target',rho_targetval,...
                        'smooth_sigma',smooth_sigmaval,'seed',seedval);
                    
                    fileName=['N_',num2str(sim.params.N),'_T_',num2str(sim.params.T),...
                        '_dt_',num2str(sim.params.dt),'_sparsity_',num2str(sim.params.sparsity),...
                        '_snr_',num2str(sim.params.snr_db),'_sigma_',num2str(sim.params.sigma_r),...
                        '_tau_',num2str(sim.params.tau),'_w_',num2str(sim.params.w),'_rhoTarget_',...
                        num2str(sim.params.rho_target),'_smoothSigma_',num2str(sim.params.smooth_sigma),...
                        '_seed_',num2str(sim.params.seed),'.xlsx'];
                    
                    %   A             : NxN ground-truth directed weights (zero diagonal)
                    %   spikes        : NxT Poisson spikes
                    %   rate_latent   : NxT latent rectified rates (a.u.)
                    %   calcium       : NxT calcium (a.u.)
                    %   F             : NxT fluorescence with noise (a.u.)
                    %   Z             : NxT z-scored fluorescence per neuron
                    %   rate_smooth   : NxT Gaussian-smoothed spike-rate (Hz)
                    %   params        : struct of all parameters
                    writematrix(sim.A,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','A');
                    writematrix(sim.spikes,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','spikes');
                    writematrix(sim.rate_latent,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','rate_latent');
                    writematrix(sim.calcium,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','calcium');
                    writematrix(sim.F,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','F');
                    writematrix(sim.Z,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','Z');
                    writematrix(sim.rate_smooth,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','rate_smooth');
                    parm=[sim.params.N,sim.params.T,sim.params.dt,sim.params.sparsity,sim.params.snr_db,...
                        sim.params.sigma_r,sim.params.tau,sim.params.w,sim.params.rho_target,...
                        sim.params.smooth_sigma,sim.params.seed];
                    writematrix(parm,fullfile(pathOut,[dataType,'_',fileName]),'Sheet','params');
             end
         end
     end
% end


function sim = simulate_var_calcium(varargin)
% SIMULATE_VAR_CALCIUM  Sparse VAR(1) -> spikes -> calcium -> fluorescence.
% Clean, stable implementation with burn-in, robust scaling, clipping, and
% explicit spike amplitude. Usage identical to your original function.

% -------------------- parse inputs --------------------
p = inputParser;
addParameter(p,'N',10);
addParameter(p,'T',2400);
addParameter(p,'dt',0.05);
addParameter(p,'sparsity',0.05);
addParameter(p,'w',0.25);
addParameter(p,'sigma_r',0.2);
addParameter(p,'rho_target',0.9);
addParameter(p,'tau',1.2);
addParameter(p,'snr_db',0);
addParameter(p,'smooth_sigma',0.2);
addParameter(p,'seed',42);
parse(p,varargin{:});
P = p.Results;
rng(P.seed);

N = P.N;
T = P.T;
dt = P.dt;

% -------------------- 1) adjacency A --------------------
A = zeros(N);
mask = rand(N) < P.sparsity;
mask(eye(N)==1) = false;           % no self edges
A(mask) = P.w .* randn(nnz(mask),1);

% stabilize: scale to desired spectral radius (use eigs when available)
if any(A(:)~=0)
    try
        er = abs(eigs(A,1,'largestabs'));
        er = er(1);
    catch
        er = max(abs(eig(A)));
    end
    if er == 0, er = 1; end
    A = (P.rho_target/er) * A;
end

% -------------------- 2) latent rate VAR(1) with burn-in --------------------
% burn_in set relative to T for robustness
burn_in = min(500, max(50, round(0.1 * T)));
T_total = T + burn_in;

rate = zeros(N, T_total);
epsR = P.sigma_r * randn(N, T_total);
rate(:,1) = epsR(:,1);   % small nonzero init

for t = 2:T_total
    rate(:,t) = A * rate(:,t-1) + epsR(:,t);
end

% discard burn-in and rectify
rate = max(rate(:, burn_in+1:end), 0);   % [N x T]

% -------------------- 3) Convert latent rates -> Poisson spikes --------------------
% Robust scaling from arbitrary "rate" units to Hz: use median of positive entries
pos_rates = rate(rate > 0);
if isempty(pos_rates)
    scale = 1;
else
    target_med_hz = 2;  % desired median firing rate (Hz)
    scale = target_med_hz / (median(pos_rates) + eps);
end
lambda = scale .* rate;    % Hz
% clip unrealistic firing rates
maxLambdaHz = 100;
lambda(lambda > maxLambdaHz) = maxLambdaHz;

% Poisson counts per bin
spikes = poissrnd(lambda * dt);   % [N x T], counts per bin

% -------------------- 4) Calcium forward model --------------------
A_spike = 1.0;                     % calcium amplitude per spike (tune)
alpha = exp(-dt / P.tau);
calcium = zeros(N, T);
for t = 2:T
    calcium(:,t) = alpha * calcium(:,t-1) + A_spike * spikes(:,t);
end

% -------------------- 5) Measurement noise -> fluorescence --------------------
F_clean = calcium;
F = zeros(N, T);
snr_lin = 10^(P.snr_db / 10);
for i = 1:N
    sig = std(F_clean(i,:));
    if sig == 0, sig = 1e-6; end
    % var_noise = var_signal / snr_lin  => noise_std = sig / sqrt(snr_lin)
    noise_std = sig / sqrt(snr_lin + eps);
    F(i,:) = F_clean(i,:) + noise_std * randn(1, T);
end

% -------------------- 6) z-score fluorescence per neuron --------------------
Z = zeros(N, T);
for i = 1:N
    mu = mean(F(i,:)); sd = std(F(i,:));
    if sd == 0, sd = 1e-6; end
    Z(i,:) = (F(i,:) - mu) / sd;
end

% -------------------- 7) Smoothed spike-rate (Gaussian kernel on spikes->Hz) ------
if P.smooth_sigma > 0
    sig_samp = max(1, round(P.smooth_sigma / dt));
    g = gausswin(6 * sig_samp + 1);
    g = g / sum(g);
    rate_smooth = zeros(N, T);
    for i = 1:N
        sHz = spikes(i,:) / dt;
        rate_smooth(i,:) = conv(sHz, g, 'same');
    end
else
    rate_smooth = spikes / dt;
end

% -------------------- package output --------------------
sim = struct();
sim.A = A;
sim.spikes = spikes;
sim.rate_latent = rate;
sim.calcium = calcium;
sim.F = F;
sim.Z = Z;
sim.rate_smooth = rate_smooth;
sim.params = P;

end