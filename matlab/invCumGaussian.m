function y = invCumGaussian(x)

% INVCUMGAUSSIAN Inverser of the cumulative Gaussian.

% IVM

y = erfinv(x*2 - 1)*2/sqrt(2);