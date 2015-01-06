clear all;
close all;
clc;

addpath('~/Dropbox/matlab/');

outfile = 'zare';
savefile = 'zare.mat';

% generate data
gendata;

loadzaredata;
opt.burnin = 10000; % main burnin phase
opt.start_it = opt.burnin + opt.initial_phase;
opt.nits = 100000 + opt.start_it; % number of MCMC sweeps

rng('default');

% use hamming ball sampling
for m = params.K : -1 : 1
    disp(['m = ' num2str(m)]);
    opt.m = m;
    out{m} = hbsample(dat, params, opt);
end

% use gibbs sampling
opt.pr_slice = 0.0;
out_gibbs = hbsample(dat, params, opt);

% plot
plotall;

save(outfile, 'out', 'out_gibbs', 'dat', 'params', 'opt');



