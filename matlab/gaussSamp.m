function y = gaussSamp(Sigma, numSamps)

% GAUSSSAMP Sample from a Gaussian with a given covariance.

% NDLUTIL

[U, V] = eig(Sigma);
dims = size(Sigma, 1);
y = randn(numSamps, dims);
y = y*diag(sqrt(diag(V)));
y = y*U';