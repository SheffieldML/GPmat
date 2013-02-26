function [dlnZ_dmu, dlnZ_dvs] = noiseGradVals(noise, mu, varsigma, y)

% NOISEGRADVALS Gradient of noise model wrt mu and varsigma.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = NOISEGRADVALS(NOISE, MU, VARSIGMA, Y)
%	computes the gradient of the given noise with respect to the input
%	mean and the input variance.
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
%	% SEEALSO NOISECREATE, NOISEPARAMINIT, NOISEGRADIENTPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence


fhandle = str2func([noise.type 'NoiseGradVals']);
[dlnZ_dmu, dlnZ_dvs] = fhandle(noise, mu, ...
                               varsigma, y);
