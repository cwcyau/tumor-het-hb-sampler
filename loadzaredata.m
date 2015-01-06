% load zare data
variantdat = importdata('~/Dropbox/external_papers/hammingball/data/tumor/17loci-variants.txt');
totaldat = importdata('~/Dropbox/external_papers/hammingball/data/tumor/17loci-total.txt');

dat.ya = variantdat.data';
dat.yd = totaldat.data';

[dat.nSamples, dat.N] = size(dat.ya);

opt.sigma_v = 1*ones(1, dat.nSamples); % variance of random walk proposal

dat.binocoeff = zeros(dat.nSamples, dat.N);
for s = 1 : dat.nSamples
	dat.binocoeff(s, :) = gammaln(dat.yd(s, :) + 1) - gammaln(dat.ya(s, :) + 1) - gammaln(dat.yd(s, :) - dat.ya(s, :) + 1);
end