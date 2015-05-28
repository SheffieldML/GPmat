function y = ngaussian(x)

% NGAUSSIAN Compute a Gaussian with mean 0 and variance 1.
% FORMAT
% DESC computes a the likelihood of a normalised Gaussian 
% distribution, i.e. with mean 0 and variance 1.
% ARG x : input value(s) for which to compute the distribution.
% RETURN y : probability of the input values under the Gaussian.
%
% SEEALSO : cumGaussian, gaussOverDiffcumGaussian
%
% COPYRIGHT : Neil D. Lawrence, 2004

% NDLUTIL

x2 = x.*x;
y = exp(-.5*x2);
y = y/sqrt(2*pi);
