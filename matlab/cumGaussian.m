function y = cumGaussian(x)

% CUMGAUSSIAN Cummulative distribution for Gaussian.

% IVM

y = 0.5*(1+myerf(sqrt(2)/2*x));
