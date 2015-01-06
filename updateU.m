function params = updateU(dat, params, opt)

for i = 1 : dat.N
    U_HB = opt.Xtable(opt.loc(params.X(i), :), :);		
    ind = randsample(1:opt.M, 1);
    params.U(i) = opt.loc(params.X(i), ind);
end	
