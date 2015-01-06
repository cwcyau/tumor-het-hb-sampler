clear all;
close all;
clc;

addpath('~/matlab/');

outfile = 'zare-conditional';
savefile = 'zare-conditional.mat';

K = 10; % number of populations
S = 2; % number of states
nConfig = S^K; % total number of column configurations

%%
%% pre-calculate quantities for hamming ball sampler
%%
D = [0 : nConfig-1]; % sequence of decimals

% convert decimal numbers to binary representations
Xtable = zeros(nConfig, K); % pre-allocate matrix of binary representations of decimals
for i = 1 : nConfig
	% for each decimal, convert to binary representation
	dec = D(i);
	Xtable(i, :) = de2bi(dec, K);
end


%%
%% simulate data
%%
%
%gendata;

loadzaredata;

rand('state', cputime);
randn('state', cputime);


%%
%% setup posterior inference using mcmc
%%

w0 = (alpha/K)*ones(1, K);
X = randsample(1:nConfig, N, 1);

gam = zeros(K, nSamples);
for s = 1 : nSamples
	gam(:, s) = 1; %gamrnd(alpha/K, 1, [K 1]);
	w(:, s) = gam(:, s)./sum(gam(:, s), 1);
end
v = log(gam);

f = f0*ones(1, nConfig);
% compute probability of each column
for i = 1 : nConfig
	Xcol = Xtable(i, :);
	pr_Xtable(i) = sum(Xcol)*log(f(i)) + (K-sum(Xcol))*log(1-f(i));
end

X_all = zeros(N, nits);
v_all = zeros(K, nSamples, nits);
w_all = zeros(K, nSamples, nits);
phi_all = zeros(nSamples, N, nits);
res_all = zeros(nSamples, nits);
loglik_all = zeros(nits, 1);

accept = zeros(nits, nSamples);
timeloop = zeros(1, nits);

for it = 1 : nits

	tic;

	if mod(it, 100) == 0
		disp(['Iteration: ' num2str(it)]);
		disp(['Accept: ' num2str(sum(accept, 1))]);
		disp(['sigma_v: ' num2str(sigma_v)]);
		tt = sum(timeloop(1:it))/it;
		disp(['Time (loop): ' num2str(tt)]);
		disp(['log-likelihood: ' num2str(loglik_all(it-1))]);
	end

	% update W conditionally
	XX = Xtable(X, :);
	loglik = zeros(1, nSamples);
	for s = 1 : nSamples
		
		v_new = v(:, s) + sigma_v(s)*randn(K, 1);
		gam_new = exp(v_new);
		w_new = gam_new./sum(gam_new);

		phi_new = 0.5*XX*w_new;
		pp = (1-e)*phi_new + e*(1-phi_new);
		loglik_new(s) = sum( binocoeff(s, :) + ya(s, :).*log(pp') + (yd(s, :)-ya(s, :)).*log(1-pp') ) + sum(loggampdf(v_new, alpha/K, 1));
		
		phi = 0.5*XX*w(:, s);
		pp = (1-e)*phi + e*(1-phi);
		loglik(s) = sum( binocoeff(s, :) + ya(s, :).*log(pp') + (yd(s, :)-ya(s, :)).*log(1-pp') ) + sum(loggampdf(v(:, s), alpha/K, 1));
		
		if log(rand) < min(0, loglik_new(s)-loglik(s))
			v(:, s) = v_new;
			w(:, s) = w_new';
			accept(it, s) = 1;
		end
		
	end	
	
	% update X
	for i = randperm(N)
		
		% compute likelihood	
		loglik_sum = pr_Xtable + obslik(ya(:, i), yd(:, i), Xtable, w, e, binocoeff(:, i), nSamples, nConfig);

		% normalise
		pX = exp( loglik_sum - logsumexp(loglik_sum, 2) );
		ind = randsample(1:nConfig, 1, 1, pX);
		
		X(i) = ind;
		
	end
	

	% compute current log-posterior
	loglik_vec = zeros(1, N);
	for i = 1 : N
		loglik_vec(i) = 0*pr_Xtable(X(i)) + obslik(ya(:, i), yd(:, i), Xtable(X(i), :), w, e, binocoeff(:, i), nSamples, 1);
	end	
	
	logprior = 0;
	for s = 1 : nSamples
		logprior = logprior + sum(loggampdf(v(:, s), alpha/K, 1));
	end
	loglik_all(it) = sum(loglik_vec) + 0*logprior; 
		

	% update proposal variance
	if it > burnin
		if mod(it, it_blk) == 0
			for s = 1 : nSamples
				accept_rate = sum(accept([it-it_blk:it-1], s), 1)/it_blk;
				if accept_rate < 0.2 & sigma_v > sigma_v_min
					sigma_v(s) = (1-sigma_v_step)*sigma_v(s);
				end
				if accept_rate > 0.4 & sigma_v < sigma_v_max
					sigma_v(s) = (1+sigma_v_step)*sigma_v(s);
				end
			end
		end
	end


	% store summary stats
	X_HB = Xtable(X, :);
	for s = 1 : nSamples		
		phi_all(s, :, it) = 0.5*X_HB*w(:, s);
%		res_all(s, it) = sqrt( mean( (X_HB*w(:, s) - Xtrue'*Wtrue(:, s)).^2 ) );
	end
	

	X_all(:, it) = X;
	w_all(:, :, it) = w;
	v_all(:, :, it) = v;

	timeloop(it) = toc;

end

% plot
plotall;

saveall;
