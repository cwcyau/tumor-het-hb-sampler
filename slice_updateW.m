function [ params, accept ] = slice_updateW(dat, params, opt)

% marginally propose w,
params_new = params;
accept = zeros(1, dat.nSamples);

for s = 1 : dat.nSamples

    if rand < opt.pr_indp_prop
        gam_new = gamrnd(params.alpha/params.K, 1, [params.K 1]);
        v_new = log(gam_new);
    else
        v_new = params.v(:, s) + opt.sigma_v(s)*randn(params.K, 1);
    end
    gam_new = exp(v_new);
    w_new = gam_new./sum(gam_new);

    % sum over Hamming Ball
    loglik_new = zeros(opt.M, dat.N);
    loglik_old = zeros(opt.M, dat.N);

    for i = 1 : dat.N

        % find X's in the Hamming Ball
        X_HB = opt.Xtable(opt.loc(params.U(i), :), :);
        %log_pr_X_HB = opt.pr_Xtable(opt.loc(params.U(i), :));
        X_HB_sum = sum(X_HB, 2);
		log_pr_X_HB = X_HB_sum.*log(params.f(i)) + (params.K-X_HB_sum).*log(1-params.f(i));

        % compute likelihood for proposed w
        loglik_new(:, i) = log_pr_X_HB + calcObslik(dat.ya(s, i), dat.yd(s, i), dat.binocoeff(s, i), X_HB, w_new, params.e);

        % compute likelihood for existing w	
        loglik_old(:, i) = log_pr_X_HB + calcObslik(dat.ya(s, i), dat.yd(s, i), dat.binocoeff(s, i), X_HB, params.w(:, s), params.e);      
        
    end
    loglik_new_sum_w = sum(logsumexp(loglik_new, 1), 2) + sum(loggampdf(v_new, params.alpha/params.K, 1));
    loglik_old_sum_w = sum(logsumexp(loglik_old, 1), 2) + sum(loggampdf(params.v(:, s), params.alpha/params.K, 1));
    
    if log(rand) < min(0, loglik_new_sum_w - loglik_old_sum_w)
        
        params.w(:, s) = w_new;
        params.gam(:, s) = gam_new;
        params.v(:, s) = v_new;
        
        accept(s) = 1;
        
    end

end
