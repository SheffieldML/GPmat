function [dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y)

% GAUSSIANNOISEGRADVALS Gradient of GAUSSIAN noise log Z with respect to input mean and variance.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = GAUSSIANNOISEGRADVALS(NOISE, MU, VARSIGMA, Y)
%	computes the gradient of the Gaussian noise with respect to the
%	input mean and the input variance.
%	 Returns:
%	  DLNZ_DMU - the gradient of log Z with respect to the input mean.
%	  DLNZ_DVS - the gradient of log Z with respect to the input
%	   variance.
%	 Arguments:
%	  NOISE - noise structure for which gradients are being computed.
%	  MU - mean input locations with respect to which gradients are
%	   being computed.
%	  VARSIGMA - variance input locations with respect to which
%	   gradients are being computed.
%	  Y - noise model output observed values associated with the given
%	   points.
%	
%
%	See also
%	% SEEALSO GAUSSIANNOISEPARAMINIT, GAUSSIANNOISEGRADIENTPARAM, NOISEGRADVALS, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);
nu = 1./(noise.sigma2+varsigma);
dlnZ_dmu = zeros(size(nu));
for i = 1:D
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i) - noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);
