function y = cumGaussian(x)

% CUMGAUSSIAN Cumulative distribution for Gaussian.
% FORMAT
% DESC computes the cumulative Gaussian distribution.
% ARG x : input value.
% RETURN p : output probability.
%
% SEEALSO : lnCumGaussian, lnDiffCumGaussian, erf
%
% COPYRIGHT : Neil D. Lawrence, 2004

% NDLUTIL

y = 0.5*(1+erf(sqrt(2)/2*x));
