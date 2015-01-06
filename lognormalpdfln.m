function y = lognormalpdfln(x, mu, sigma)

y = - log(x*sigma*sqrt(2*pi)) - 0.5* ((log(x) - mu)/sigma).^2 ;
