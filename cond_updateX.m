function params = cond_updateX(dat, params, opt)
%
% conditionally update X
%
for i = randperm(dat.N)
    
    % compute likelihood
    %loglik_sum = opt.pr_Xtable;
    X_sum = sum(opt.Xtable(params.X(i), :));
	loglik_sum = X_sum.*log(params.f(i)) + (params.K-X_sum).*log(1-params.f(i));

    for s = 1 : dat.nSamples
        loglik_sum = loglik_sum + calcObslik(dat.ya(s, i), dat.yd(s, i), dat.binocoeff(s, i), opt.Xtable, params.w(:, s), params.e);
    end
    % normalise
    pX = exp( loglik_sum - logsumexp(loglik_sum, 1) );
    ind = randsample(1:opt.nConfig, 1, 1, pX);
    
    params.X(i) = ind;
    
end