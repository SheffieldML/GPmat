function f = lnDiffCumGaussian(u, uprime)

% LNDIFFCUMGAUSSIAN Computes the log of the difference between two cumulative Gaussians.

% NDLUTIL

% f = log(\phi(u) - phi(uprime))


f = log(gaussOverDiffCumGaussian(u, uprime, 1)+1e-300) + .5*u.*u + .5*log(2*pi);