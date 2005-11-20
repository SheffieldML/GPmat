function y = invCumGaussian(x)

% INVCUMGAUSSIAN Inverser of the cumulative Gaussian.

% NDLUTIL

y = erfinv(x*2 - 1)*2/sqrt(2);