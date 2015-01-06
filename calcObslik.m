function loglik = calcObslik(ya, yd, binocoeff, X_HB, w, e)

% compute variant allele frequency
phi = 0.5*X_HB*w;
pp = (1-e)*phi + e*(1-phi);

% compute log-likelihood
loglik = binocoeff + ya.*log(pp) + (yd-ya).*log(1-pp);
