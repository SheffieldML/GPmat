function kl = kldivGaussian(mean1, cov1, mean2, cov2)

% KLDIVGAUSSIAN Give the KL divergence between two Gaussians.
% FORMAT
% DESC returns the Kullback-Leibler divergence between two
% Gaussians with given means and covariances.
% ARG mean1 : mean of the first Gaussian.
% ARG cov1 : covariance of the first Gaussian.
% ARG mean2 : mean of the second Gaussian.
% ARG cov2 : covariance of the second Gaussian.
%
% SEEALSO : logdet, pdinv
%
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL

[invCov2, U] = pdinv(cov2);
logDet2 = logdet(cov2, U);
logDet1 = logdet(cov1);
N = size(cov1, 1);
meanDiff = mean1 - mean2;
if size(meanDiff, 1) == 1
  meanDiff = meanDiff';
end
kl = -0.5*(logDet1 - logDet2 - trace(cov1*invCov2) ...
           + N - meanDiff'*invCov2*meanDiff);

