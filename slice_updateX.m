function params = slice_updateX(dat, params, opt)
%
% update X marginalising over the Hamming Ball U
%
for i = randperm(dat.N)
    
    % find X's in the Hamming Ball
    X_HB = opt.Xtable(opt.loc(params.U(i), :), :);
    %log_pr_X_HB = opt.pr_Xtable(opt.loc(params.U(i), :));
    X_HB_sum = sum(X_HB, 2);
	log_pr_X_HB = X_HB_sum.*log(params.f(i)) + (params.K-X_HB_sum).*log(1-params.f(i));
    
    % compute likelihood
    loglik_sum = log_pr_X_HB;
    for s = 1 : dat.nSamples
        loglik_sum = loglik_sum + calcObslik(dat.ya(s, i), dat.yd(s, i), dat.binocoeff(s, i), X_HB, params.w(:, s), params.e);
    end
    
    % normalise
    pX = exp( loglik_sum - logsumexp(loglik_sum, 1) );
    ind = randsample(1:opt.M, 1, 1, pX);
    
    params.X(i) = opt.loc(params.U(i), ind);
    
end
