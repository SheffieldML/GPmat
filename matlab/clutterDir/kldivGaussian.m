function kl = kldivGaussian(mean1, cov1, mean2, cov2)

% KLDIVGAUSSIAN Give the KL divergence between two Gaussians.

% NDLUTILS

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

