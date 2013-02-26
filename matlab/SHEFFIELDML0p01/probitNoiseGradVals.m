function [dlnZ_dmu, dlnZ_dvs] = probitNoiseGradVals(noise, mu, varsigma, y)

% PROBITNOISEGRADVALS Gradient of PROBIT noise log Z with respect to input mean and variance.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = PROBITNOISEGRADVALS(NOISE, MU, VARSIGMA, Y)
%	computes the gradient of the probit based classification noise with
%	respect to the input mean and the input variance.
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
%	% SEEALSO PROBITNOISEPARAMINIT, PROBITNOISEGRADIENTPARAM, NOISEGRADVALS, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(mu, 2);
c = y./sqrt(noise.sigma2+varsigma);
u = zeros(size(c));
for i = 1:D
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
end
dlnZ_dmu = c.*gradLogCumGaussian(u);
dlnZ_dvs = -.5*c.*u.*dlnZ_dmu;
