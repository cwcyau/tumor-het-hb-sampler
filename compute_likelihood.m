function loglik_all = compute_likelihood(dat, params, opt);

loglik_vec = zeros(dat.nSamples, dat.N);
for s = 1 : dat.nSamples
    for i = 1 : dat.N
        X_HB = opt.Xtable(params.X(i), :);
        loglik_vec(s, i) = calcObslik(dat.ya(s, i), dat.yd(s, i), dat.binocoeff(s, i), X_HB, params.w(:, s), params.e);
    end	
end

logprior = 0;
for s = 1 : dat.nSamples
    logprior = logprior + sum(loggampdf(params.v(:, s), params.alpha/params.K, 1));
end
loglik_all = sum(loglik_vec(:)) + 0*logprior; 
