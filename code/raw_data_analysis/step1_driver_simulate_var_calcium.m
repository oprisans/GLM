% step1_driver_simulate_var_calcium.m  (fixed)
% Parallelized over seedList. Writes ground-truth .xlsx files to data/groundTruth.
%
% Usage examples:
%   seedList = [42 101 202 303];
%   parpool('local',6);
%   run('driver_simulate_var_calcium.m');

%% ----- project paths & defaults (caller may override these variables in base) -----
scriptFull = mfilename('fullpath');
if isempty(scriptFull)
    codeDir = pwd;
else
    codeDir = fileparts(scriptFull);
end
projectRoot = fileparts(codeDir);    % .../simulations2

% output folder
dataPath    = fullfile(projectRoot,'data');
outGTdir    = fullfile(dataPath,'groundTruth');
if ~isfolder(outGTdir), mkdir(outGTdir); end

% defaults (caller may override in base workspace)
if ~exist('Nvals','var'),        Nvals = 10:10:30; end
if ~exist('Tval','var'),         Tval = 2400; end
if ~exist('dtval','var'),        dtval = 0.05; end
if ~exist('sparsityList','var'), sparsityList = 0.05:0.05:0.30; end
if ~exist('snr_vals','var'),     snr_vals = 0:5:10; end
if ~exist('wval','var'),         wval = 0.25; end
if ~exist('sigma_rval','var'),   sigma_rval = 0.2; end
if ~exist('rho_targetval','var'),rho_targetval = 0.9; end
if ~exist('tauval','var'),       tauval = 1.2; end
if ~exist('smooth_sigmaval','var'), smooth_sigmaval = 0.2; end

% seed handling (caller may set seedList or single seedval)
if ~exist('seedList','var') || isempty(seedList)
    if exist('seedval','var') && ~isempty(seedval)
        seedList = seedval(:).';
    else
        seedList = 42;
    end
end

% ensure pool workers see the code folder
pool = gcp('nocreate');
if ~isempty(pool)
    pctRunOnAll addpath(codeDir);
end

% parallelize over seeds: each worker sets its RNG once, then produces multiple datasets
parfor si = 1:numel(seedList)
    s = seedList(si);

    % ---- deterministic per-worker stream ----
    % use explicit RandStream object and set it as the worker global stream.
    % mt19937ar is default and portable. This is set once per worker here.
    rs = RandStream('mt19937ar','Seed', s);
    RandStream.setGlobalStream(rs);

    % now draw reproducible random numbers on this worker.
    for Ni = 1:numel(Nvals)
        N = Nvals(Ni);
        for spi = 1:numel(sparsityList)
            sparsity = sparsityList(spi);
            for svi = 1:numel(snr_vals)
                snr_db = snr_vals(svi);

                try
                    % call the helper. DO NOT reseed inside the helper.
                    sim = simulate_var_calcium('N',N,'T',Tval,'dt',dtval,'sparsity',sparsity,...
                        'snr_db',snr_db,'sigma_r',sigma_rval,'tau',tauval,'w',wval,'rho_target',rho_targetval,...
                        'smooth_sigma',smooth_sigmaval,'seed',s);

                    fileName = sprintf('groundTruth_N_%d_T_%d_dt_%.5g_sparsity_%.3g_snr_%g_sigma_%.3g_tau_%.3g_w_%.3g_rhoTarget_%.3g_smoothSigma_%.3g_seed_%d.xlsx', ...
                        sim.params.N, sim.params.T, sim.params.dt, sim.params.sparsity, sim.params.snr_db, ...
                        sim.params.sigma_r, sim.params.tau, sim.params.w, sim.params.rho_target, sim.params.smooth_sigma, sim.params.seed);

                    outFile = fullfile(outGTdir, fileName);

                    % Write sheets
                    writematrix(sim.A, outFile, 'Sheet', 'A');
                    writematrix(sim.spikes, outFile, 'Sheet', 'spikes');
                    writematrix(sim.rate_latent, outFile, 'Sheet', 'rate_latent');
                    writematrix(sim.calcium, outFile, 'Sheet', 'calcium');
                    writematrix(sim.F, outFile, 'Sheet', 'F');
                    writematrix(sim.Z, outFile, 'Sheet', 'Z');
                    writematrix(sim.rate_smooth, outFile, 'Sheet', 'rate_smooth');
                    parm = [sim.params.N,sim.params.T,sim.params.dt,sim.params.sparsity,sim.params.snr_db, ...
                        sim.params.sigma_r,sim.params.tau,sim.params.w,sim.params.rho_target, ...
                        sim.params.smooth_sigma,sim.params.seed];
                    writematrix(parm, outFile, 'Sheet', 'params');

                catch ME
                    warning('seed %d N=%d sparsity=%.3g snr=%g failed: %s', s, N, sparsity, snr_db, ME.message);
                end

            end
        end
    end
end

%% ----------------- local helper: simulate_var_calcium -----------------
function sim = simulate_var_calcium(varargin)
% SIMULATE_VAR_CALCIUM  Sparse VAR(1) -> spikes -> calcium -> fluorescence.
% Accepts parameter-value pairs identical to calls above.
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

% NOTE: DO NOT reseed here. caller (parfor worker) sets RNG once.
% rng(P.seed);   <-- intentionally omitted

N = P.N; T = P.T; dt = P.dt;

% adjacency
A = zeros(N);
mask = rand(N) < P.sparsity;
mask(eye(N)==1) = false;
A(mask) = P.w .* randn(nnz(mask),1);

% stabilize spectral radius
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

% VAR(1) latent rate with burn-in
burn_in = min(500, max(50, round(0.1 * T)));
T_total = T + burn_in;
rate = zeros(N, T_total);
epsR = P.sigma_r * randn(N, T_total);
rate(:,1) = epsR(:,1);
for t = 2:T_total
    rate(:,t) = A * rate(:,t-1) + epsR(:,t);
end
rate = max(rate(:, burn_in+1:end), 0);

% scale to Hz and Poisson spikes
pos_rates = rate(rate > 0);
if isempty(pos_rates)
    target_med_hz = 2;
    scale = 1;
else
    target_med_hz = 2;
    scale = target_med_hz / (median(pos_rates) + eps);
end
lambda = scale .* rate;
maxLambdaHz = 100;
lambda(lambda > maxLambdaHz) = maxLambdaHz;
spikes = poissrnd(lambda * dt);

% calcium dynamics
A_spike = 1.0;
alpha = exp(-dt / P.tau);
calcium = zeros(N, T);
for t = 2:T
    calcium(:,t) = alpha * calcium(:,t-1) + A_spike * spikes(:,t);
end

% fluorescence with additive Gaussian noise (SNR in dB)
F_clean = calcium;
F = zeros(N, T);
snr_lin = 10^(P.snr_db / 10);
for i = 1:N
    sig = std(F_clean(i,:));
    if sig == 0, sig = 1e-6; end
    noise_std = sig / sqrt(snr_lin + eps);
    F(i,:) = F_clean(i,:) + noise_std * randn(1, T);
end

% z-score per neuron
Z = zeros(N, T);
for i = 1:N
    mu = mean(F(i,:)); sd = std(F(i,:));
    if sd == 0, sd = 1e-6; end
    Z(i,:) = (F(i,:) - mu) / sd;
end

% Gaussian smoothing of spike-rate (Hz)
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

% package output
sim.A = A;
sim.spikes = spikes;
sim.rate_latent = rate;
sim.calcium = calcium;
sim.F = F;
sim.Z = Z;
sim.rate_smooth = rate_smooth;
sim.params = P;

end % end of simulate_var_calcium