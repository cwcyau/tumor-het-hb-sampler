clear all;
close all;
clc;

addpath('~/Dropbox/matlab/');

outfile = 'longrun-10';
savefile = 'longrun-10.mat';

% generate data
gendata;

loadzaredata

rng('default');

% use hamming ball sampling
opt.m = 4;
out = hbsample(dat, params, opt);

% use gibbs sampling
%opt.pr_slice = 0.0;
%out_gibbs = hbsample(dat, params, opt);

% plot
plotall;

%save(outfile, 'out', 'out_gibbs', 'dat', 'params', 'opt');

