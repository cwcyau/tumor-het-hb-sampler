function y = gampdfln(x, k, theta)
%
% p = x^(k-1) exp(-x/theta) / theta^k Gamma(k)
%
% k - shape
% theta - scale
%
y = (k-1)*log(x) - x./theta - k*log(theta) - gammaln(k);
