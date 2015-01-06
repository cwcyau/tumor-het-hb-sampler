function out = hbsample(dat, params, opt)

%%
%% pre-calculate quantities for hamming ball sampler
%%
D = [ 0 : opt.nConfig-1 ]; % sequence of decimals

% convert decimal numbers to binary representations
opt.Xtable = zeros(opt.nConfig, params.K); % pre-allocate matrix of binary representations of decimals
for i = 1 : opt.nConfig
	% for each decimal, convert to binary representation
	dec = D(i);
	opt.Xtable(i, :) = de2bi(dec, params.K);
end

% compute hamming distance between all numbers
dist = zeros(opt.nConfig, opt.nConfig); % hamming distance matrix 
for i = 1 : opt.nConfig
	tmp = bsxfun(@xor, opt.Xtable(i, :), opt.Xtable);
	dist(:, i) = sum(tmp, 2);
end

% locations of vectors satisfying hamming distance criterion
opt.M = 0;
for j = 0 : opt.m
	opt.M = opt.M + nchoosek(params.K, j);
end
opt.loc = zeros(opt.nConfig, opt.M);
for i = 1 : opt.nConfig
	opt.loc(i, :) = find( dist(:, i) <= opt.m );
end

%%
%% setup posterior inference using mcmc
%%
params.w0 = (params.alpha/params.K)*ones(1, params.K);
params.X = randsample(1:opt.nConfig, dat.N, 1);
params.U = params.X;

params.gam = zeros(params.K, dat.nSamples);
for s = 1 : dat.nSamples
	params.gam(:, s) = gamrnd(params.alpha/params.K, 1, [params.K 1]);
	params.w(:, s) = params.gam(:, s)./sum(params.gam(:, s), 1);
end
params.v = log(params.gam);
params.f = params.f0*ones(1, dat.N);

% compute probability of each column
opt.pr_Xtable = zeros(opt.nConfig, 1);
for i = 1 : opt.nConfig
	Xcol = opt.Xtable(i, :);
	opt.pr_Xtable(i) = sum(Xcol)*log(params.f0) + (params.K-sum(Xcol))*log(1-params.f0);
end

out.X_all = zeros(dat.N, opt.nits, 'uint32');
out.v_all = zeros(params.K, dat.nSamples, opt.nits, 'single');
out.w_all = zeros(params.K, dat.nSamples, opt.nits);
out.f_all = zeros(dat.N, opt.nits, 'single');
out.phi_all = zeros(dat.nSamples, dat.N, opt.nits, 'single');
out.res_all = zeros(dat.nSamples, opt.nits, 'single');
out.loglik_all = zeros(opt.nits, dat.nSamples);
out.accept = zeros(opt.nits, dat.nSamples, 'uint8');
out.timeloop = zeros(1, opt.nits);

for it = 1 : opt.nits

	tic;

	if mod(it, opt.n_time_blk) == 0
		disp(['Iteration: ' num2str(it)]);
        if it > opt.initial_phase
            disp(['Accept: ' num2str(sum(out.accept(it-opt.it_blk:it, :), 1)/opt.it_blk)]);
        else
            disp(['Accept: ' num2str(sum(out.accept, 1))]);
        end
		disp(['sigma_v: ' num2str(opt.sigma_v)]);
		tt = sum(out.timeloop(1:it))/it;
		disp(['Time (loop): ' num2str(tt)]);
		disp(['log-likelihood: ' num2str(out.loglik_all(it-1))]);
        for s = 1 : dat.nSamples
            disp(['w(' num2str(s) ') = ' num2str(params.w(:, s)')]);
        end
		disp(['f = ' num2str(params.f)]);
    end

    if rand < opt.pr_slice

        % slice update W
        [ params, out.accept(it, :) ] = slice_updateW(dat, params, opt);
        
        % slice update X
        params = slice_updateX(dat, params, opt);
        
    else

        % update W conditionally
        [ params, out.accept(it, :) ] = cond_updateW(dat, params, opt);
        
        % update X conditionally
        params = cond_updateX(dat, params, opt);

    end

	% move U
    params = updateU(dat, params, opt);
    
    % update frequencies
    Xsum = sum(opt.Xtable(params.X, :), 2);
    for i = 1 : dat.N
    	params.f(i) = betarnd(params.f_alpha + Xsum(i), params.f_beta + (params.K-Xsum(i)), [1 1]);
    end

	% update proposal variance during tuning phase
	if it > opt.initial_phase & it < opt.start_it
		if mod(it, opt.it_blk) == 0     
            opt = update_proposal_variance(out.accept([it-opt.it_blk:it-1], :), dat, params, opt);
		end
    end
	
	% compute overall likelihood
    out.loglik_all(it, :) = compute_likelihood(dat, params, opt);	
		
	% store output
	out.X_all(:, it) = params.X;
	out.w_all(:, :, it) = params.w;
	out.v_all(:, :, it) = params.v;
	out.f_all(:, it) = params.f;
	
	% compute residuals
	X_mat = opt.Xtable(params.X, :);
	for s = 1 : dat.nSamples		
		out.phi_all(s, :, it) = 0.5*X_mat*params.w(:, s);
%		out.res_all(s, it) = sqrt( mean( (X_mat*w(:, s) - Xtrue'*Wtrue(:, s)).^2 ) );
    end

	out.timeloop(it) = toc;

end

out.params = params;
out.opt = opt;
