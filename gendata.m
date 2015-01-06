
%%
%% simulate data
%%

rng('default');
rng(1);

params.K = 10; % number of populations
opt.S = 2; % number of states
opt.nConfig = opt.S^params.K; % total number of column configurations
opt.m = 2; % max. hamming distances


%% MCMC options
%%
opt.pr_indp_prop = 0.1; % independent proposal
opt.pr_slice = 1; % probability of slice sampling proposal (0 - Conditional Gibbs Update, 1 - Hamming Ball Update)
opt.initial_phase = 100; % initial tuning phase
opt.burnin = 1000; % main burnin phase
opt.start_it = opt.burnin + opt.initial_phase;
opt.nits = 10000 + opt.start_it; % number of MCMC sweeps
opt.sigma_v_step = 0.1; % step-size of random walk adjustment
opt.sigma_v_min = 0.01; % minimum allowed variance of random walk proposal
opt.sigma_v_max = 10; % maximum allowed variance of random walk proposal
opt.it_blk = 50; % number of MCMC sweeps before updating variance proposal
opt.accept_lower_threshold = 0.1; % low acceptance rate threshold (acceptance rate is adjusted to be between the lower and upper values)
opt.accept_upper_threshold = 0.3; % upper acceptance rate threshold (acceptance rate is adjusted to be between the lower and upper values)
opt.n_time_blk = 100; % number of iterations before printing timing information

%% hyperpararameters
params.alpha = 1; % dirichlet prior parameter
params.f0 = 0.5; % base frequency of ones
params.f_alpha = 0.5;
params.f_beta = 0.5;
params.e = 0.001; % read error probability
		 
%% simulation parameters
dat.nSamples = 1; % number of samples
dat.depthCoverage = 800; 
opt.sigma_v = 1*ones(1, dat.nSamples); % variance of random walk proposal

reps = 10;
Xtrue = [ repmat([ 1 1 1 ]', [1 reps]) repmat([ 1 1 0 ]', [1 reps]) repmat([ 1 0 0 ]' , [1 reps]) ] ;

[K0, dat.N] = size(Xtrue);		  
		  
Wtrue(:, 1) = [ 0.3 0.3 0.4 ]';
Wtrue(:, 2) = [ 0.55 0.15 0.3 ]';
Wtrue(:, 3) = [ 0.0 0.0 1.0 ]';

dat.yd = zeros(dat.nSamples, dat.N);
dat.ya = zeros(dat.nSamples, dat.N);
for s = 1 : dat.nSamples
	dat.yd(s, :) = dat.depthCoverage; %poissrnd(dat.depthCoverage, 1, dat.N);
	dat.ya(s, :) = binornd(dat.yd(s, :)', 0.5*Xtrue'*Wtrue(:, s))';
end

dat.binocoeff = zeros(dat.nSamples, dat.N);
for s = 1 : dat.nSamples
	dat.binocoeff(s, :) = gammaln(dat.yd(s, :) + 1) - gammaln(dat.ya(s, :) + 1) - gammaln(dat.yd(s, :) - dat.ya(s, :) + 1);
end
