function opt = update_proposal_variance(accept, dat, params, opt)

for s = 1 : dat.nSamples
    accept_rate = sum(accept(:, s), 1)/opt.it_blk;
    if accept_rate < opt.accept_lower_threshold & opt.sigma_v(s) > opt.sigma_v_min
        opt.sigma_v(s) = (1-opt.sigma_v_step)*opt.sigma_v(s);
    end
    if accept_rate > opt.accept_upper_threshold & opt.sigma_v(s) < opt.sigma_v_max
        opt.sigma_v(s) = (1+opt.sigma_v_step)*opt.sigma_v(s);
    end
end