function f = lnDiffCumGaussian(u, uprime)

% LNDIFFCUMGAUSSIAN Log of the difference between two cumulative Gaussians.
% FORMAT
% DESC computes the logarithm of the difference between two
% cumulative Gaussian distributions.
% ARG u1 : the argument of the first (positive) cumulative
% Gaussian.
% ARG u2 : the argument of the second (negative) cumulative
% Gaussian.
% RETURN f : where f = log(cumGaussian(u1) - cumGaussian(u2)).
%
% SEEALSO : cumGaussian, gaussOverDiffCumGaussian, lnCumGaussian
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% NDLUTIL

f = log(gaussOverDiffCumGaussian(u, uprime, 1)+1e-300) ...
    + .5*u.*u + .5*log(2*pi);
f=-f;
