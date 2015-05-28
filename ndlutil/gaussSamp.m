function y = gaussSamp(Sigma, numSamps)

% GAUSSSAMP Sample from a Gaussian with a given covariance.
% FORMAT 
% DESC samples a given number of samples from a Gaussian with a
% given covariance matrix.
% ARG Sigma : the covariance of the Gaussian to sample from.
% ARG numSamps : the number of samples to take from Gaussian.
% RETURN y : the samples from the Gaussian
%
% SEEALSO : randn, eig
%
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL

[U, V] = eig(Sigma);
dims = size(Sigma, 1);
y = randn(numSamps, dims);
y = y*diag(sqrt(diag(V)));
y = y*U';
