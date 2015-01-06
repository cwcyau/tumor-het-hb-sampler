clear all;
close all;
clc;

addpath('~/matlab/');

K = 5; % number of populations
S = 2; % number of states
N = S^K; % total number of column configurations
m = 1; % max. hamming distances


%%
%% pre-calculate quantities for hamming ball sampler
%%

D = [0 : N-1]; % sequence of decimals

% convert decimal numbers to binary representations
X = zeros(N, K); % pre-allocate matrix of binary representations of decimals
for i = 1 : N
	% for each decimal, convert to binary representation
	dec = D(i);
	X(i, :) = bitget(dec, 1:K);
end

% compute hamming distance between all numbers
dist = zeros(N, N); % hamming distance matrix 
for i = 1 : N
	tmp = bsxfun(@xor, X(i, :), X);
	dist(:, i) = sum(tmp, 2);
end

% locations of vectors satisfying hamming distance criterion
M = 0;
for j = 0 : m
	M = M + nchoosek(K, j);
end
loc = zeros(N, M);
for i = 1 : N
	loc(i, :) = find( dist(:, i) <= m );
end


%%
%% simulate data
%%

% define parameters
T = 500; % length of seq
rho = 0.05; % switching prob
sigma = 0.1; % obs. noise std. dev.
sigma2 = sigma^2;
h = 2.5; % linear term
w = zeros(1, K);
w(1:3) = [ 0.5 0.3 0.2 ]; % population weights

% generate transition matrix
Pi = zeros(N, N);
for i = 1 : N
	for j = 1 : N
		if i==j
			Pi(i, j) = 1-rho;
		else
			Pi(i, j) = rho/(N-1);
		end
	end
end

% generate simulated latent sequence
x = ones(T, 1);
%x(1) = randsample(1:N, 1);
%for t = 2 : T
%	x(t) = randsample(1:N, 1, 1, Pi(x(t-1), :));
%end
x(100:200) = 2;
x(400:500) = 7;
x(1:50) = 3;
x(300:350) = 5;

% generate simulated observed sequence
y = sum(repmat(w, [T 1]).*X(x, :), 2)*h + sigma*randn(T, 1);
x0 = x; % store true latent sequence


%%
%% setup posterior inference using mcmc
%%

nits = 1000; % number of MCMC sweeps

% preallocate storage matrices for mcmc output
prevalence_all = zeros(nits, T);
obslike_X = zeros(T, N);
obslike_X_new = zeros(T, N);
w_all = zeros(nits, K);
x_all = zeros(nits, T);
y_all = zeros(nits, T);
rho_all = zeros(nits, 1);
h_all = zeros(nits, 1);
loglik_all = zeros(nits, 1);
timeTaken = zeros(nits, 1);

% setup priors
rho = 0.01; % initialise switching prob
rho_a = 0.5; % rho ~ Beta(rho_a, rho_b) 
rho_b = 0.5; %

% setup transition matrix
Pi = zeros(N, N);
for i = 1 : N
	for j = 1 : N
		if i==j
			Pi(i, j) = 1-rho;
		else
			Pi(i, j) = rho/(N-1);
		end
	end
end

% h ~ Normal(h0, h0_sigma^2)
h0 = 0; 
h0_sigma = 5;
h0_sigma2 = h0_sigma^2;

% weights defined using normalised gamma representaton, w_k = v_k/sum_j v_j
alpha = 1/K; % w ~ Dirichlet(alpha/K, ..., alpha/K)
v = gamrnd(alpha, 1, [1 K]);
w = v/sum(v);
v_step = 1; % step size for normal kernel proposal density

% initialisation
% randomly initialise X, U
x = randsample(1:N, T, 1);
u = x;

% calculate log-likelihood for initial params
Xsum = sum(X.*repmat(w, [N 1]), 2);
for i = 1 : N
	obslike_X(:, i) = normpdf(y, h*Xsum(i), sigma);
end
[phi, phi_sum] = forwardhb(Pi, obslike_X, u, loc, M, T);
loglik_old = sum(log(phi_sum));

% MCMC
for it = 1 : nits

	tic;

	if mod(it, 100) == 0
		disp(it);
		disp(w);
	end

	% sample uniformly from within Hamming ball of x to generate u
	rn = randsample(1:M, T, 1);
	I = sub2ind([N M], x, rn);
	u = loc(I);
	
	% sample w

	% pre-compute log-likelihood (new)

	% generate proposal using log-normal proposal
	ln_v = log(v);
	ln_v_new = ln_v + v_step*randn(1, K); 
	v_new = exp(ln_v_new);
	w_new = v_new./sum(v_new);	
	
	[w_new, I] = sort(w_new, 'descend');
	v_new = v_new(I);
	
	Xsum = sum(X.*repmat(w_new, [N 1]), 2);
	for i = 1 : N
		obslike_X_new(:, i) = normpdf(y, h*Xsum(i), sigma);
	end
		
	[phi_new, phi_sum_new] = forwardhb(Pi, obslike_X_new, u, loc, M, T);
	loglik_new = sum(log(phi_sum_new));

	% calculate log proposal probability using log-normal distribution
	log_prop_new = sum( lognormalpdfln(v_new, v, v_step) );
	log_prop_old = sum( lognormalpdfln(v_new, v, v_step) );

	% calculate log prior probability using log-normal distribution
	log_prior_new = sum(gampdfln(v_new, alpha, 1));
	log_prior_old = sum(gampdfln(v, alpha, 1));

	% calculate acceptance prob
	if log(rand) < min(loglik_new+log_prior_new-loglik_old-log_prior_old+log_prop_old-log_prop_new, 0)
		w = w_new;
		v = v_new;
		obslike_X = obslike_X_new;
		phi = phi_new;
		loglik_old = loglik_new;
	end
	
	% do forward filtering - backward sampling to sample X
	[phi, phi_sum] = forwardhb(Pi, obslike_X, u, loc, M, T);
	x = backwardsamphb(Pi, phi, u, loc, M, T);

	% sample rho
	rho = betarnd(rho_a + sum(sum(x(:, 2:T)~=x(1:T-1), 1), 2), rho_b + K*(T-1));
	Pi = zeros(N, N);
	for i = 1 : N
		for j = 1 : N
			if i==j
				Pi(i, j) = 1-rho;
			else
				Pi(i, j) = rho/(N-1);
			end
		end
	end

%	% sample h
%	X_sum = sum(repmat(w, [T 1]).*X(x, :), 2);
%	h_sigma2 = 1/(sum(X_sum.^2)/sigma2 + 1/h0_sigma2);
%	h_sigma = sqrt(h_sigma2);
%	h_mu = h_sigma2*( h0/h0_sigma2 + sum(X_sum.*y)/sigma2 );
%	h = h_sigma*randn + h_mu;

	% store
	rho_all(it) = rho;
	h_all(it) = h;
	[w_all(it, :), I] = sort(w, 'descend');
	x_all(it, :) = x;
	ypred = h*sum(repmat(w, [T 1]).*X(x, :), 2);
	y_all(it, :) = ypred;
	prevalence_all(it, :) = sum(repmat(w, [T 1]).*X(x, :), 2);
	
	loglik_all(it) = sum(lognormpdf(y, ypred, sigma));

	timeTaken(it) = toc;

end



% plot data
figure(1); clf;

set(gcf, 'Renderer', 'Painters');

subplot(5, 1, 1);
hold on;
plot(y, 'k.');
plot(mean(y_all, 1), 'r.');
xlim([0 T]);

subplot(5, 1, 2);
plot(sum(X(x0, :), 2), 'k.');
ylim([-0.5 K+0.5]);

subplot(5, 1, 3);
imagesc(X(x0, :)');
xlim([0 T]);

subplot(5, 1, 4);
imagesc(X(x, :)');
xlim([0 T]);

subplot(5, 1, 5);
plot(mean(prevalence_all, 1), 'k.');
xlim([0 T]);

print -dpsc2 -r300 hammingball.ps;

figure(2); clf;

set(gcf, 'Renderer', 'Painters');

subplot(2, 2, 1);
hold on;
plot(rho_all, 'k.');
xlim([0 nits]);

subplot(2, 2, 2);
hold on;
plot(h_all, 'k.');
xlim([0 nits]);

subplot(2, 2, 3);
hold on;
plot(loglik_all, 'k.');
xlim([0 nits]);

subplot(2, 2, 4);
plot( mean(w_all, 1), 'rx');
xlim([0 K+1]);
ylim([0 1]);

print -dpsc2 -r300 -append hammingball.ps;


gzip hammingball.ps;
delete hammingball.ps;
