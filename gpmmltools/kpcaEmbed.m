function [X, sigma2] = kpcaEmbed(Y, dims)

% KPCAEMBED Embed data set with kernel PCA.

% MLTOOLS


if any(any(isnan(Y)))
  error('When missing data is present Kernel PCA cannot be used to initialise')
end

K = kernCompute(kern, Y);
[u, v] = eigs(K, dims);
X = u*sqrt(v);
sigma2 = -1;
