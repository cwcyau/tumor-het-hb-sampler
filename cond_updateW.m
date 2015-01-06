function [ params, accept ] = cond_updateW(dat, params, opt)

XX = opt.Xtable(params.X, :);
loglik = zeros(1, dat.nSamples);
loglik_new = zeros(1, dat.nSamples);

accept = zeros(1, dat.nSamples);

for s = 1 : dat.nSamples

    v_new = params.v(:, s) + opt.sigma_v(s)*randn(params.K, 1);
    gam_new = exp(v_new);
    w_new = gam_new./sum(gam_new);

    phi_new = 0.5*XX*w_new;
    pp = (1-params.e)*phi_new + params.e*(1-phi_new);
    loglik_new(s) = sum( dat.binocoeff(s, :) + dat.ya(s, :).*log(pp') + (dat.yd(s, :)-dat.ya(s, :)).*log(1-pp') ) + sum(loggampdf(v_new, params.alpha/params.K, 1));

    phi = 0.5*XX*params.w(:, s);
    pp = (1-params.e)*phi + params.e*(1-phi);
    loglik(s) = sum( dat.binocoeff(s, :) + dat.ya(s, :).*log(pp') + (dat.yd(s, :)-dat.ya(s, :)).*log(1-pp') ) + sum(loggampdf(params.v(:, s), params.alpha/params.K, 1));

    if log(rand) < min(0, loglik_new(s)-loglik(s))
        params.v(:, s) = v_new;
        params.gam(:, s) = gam_new;
        params.w(:, s) = w_new';
        accept(s) = 1;
    end

end	