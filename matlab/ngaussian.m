function y = ngaussian(x)

% NGAUSSIAN Compute a Gaussian with mean 0 and variance 1.

% IVM

x2 = x.*x;
y = exp(-.5*x2);
y = y/sqrt(2*pi);
