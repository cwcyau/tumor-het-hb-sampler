function p = loggampdf(x, k, theta)
%
% logarithm of the log-gamma pdf
%
p = k.*x - exp(x)./theta - k*log(theta) - gammaln(k);
