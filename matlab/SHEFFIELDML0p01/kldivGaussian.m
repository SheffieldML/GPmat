function kl = kldivGaussian(mean1, cov1, mean2, cov2)

% KLDIVGAUSSIAN Give the KL divergence between two Gaussians.
%
%	Description:
%
%	KLDIVGAUSSIAN(MEAN1, COV1, MEAN2, COV2) returns the Kullback-Leibler
%	divergence between two Gaussians with given means and covariances.
%	 Arguments:
%	  MEAN1 - mean of the first Gaussian.
%	  COV1 - covariance of the first Gaussian.
%	  MEAN2 - mean of the second Gaussian.
%	  COV2 - covariance of the second Gaussian.
%	
%
%	See also
%	LOGDET, PDINV


%	Copyright (c) 2005 Neil D. Lawrence


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

