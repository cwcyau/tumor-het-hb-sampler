clear all;
close all;
clc;

addpath('~/Dropbox/matlab/');

outfile = 'longrun-5';
savefile = 'longrun-5.mat';

% generate data
gendata;

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

