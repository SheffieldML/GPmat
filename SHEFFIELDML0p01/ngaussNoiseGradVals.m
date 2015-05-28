function [dlnZ_dmu, dlnZ_dvs] = ngaussNoiseGradVals(noise, mu, varsigma, y)

% NGAUSSNOISEGRADVALS Gradient of NGAUSS noise log Z with respect to input mean and variance.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = NGAUSSNOISEGRADVALS(NOISE, MU, VARSIGMA, Y)
%	computes the gradient of the noiseless Gaussian noise with respect
%	to the input mean and the input variance.
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
%	% SEEALSO NGAUSSNOISEPARAMINIT, NGAUSSNOISEGRADIENTPARAM, NOISEGRADVALS, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



[dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y);