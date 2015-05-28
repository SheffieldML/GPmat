function [dlnZ_dmu, dlnZ_dvs] = cmpndNoiseGradVals(noise, mu, varsigma, y)

% CMPNDNOISEGRADVALS Gradient of CMPND noise log Z with respect to input mean and variance.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = CMPNDNOISEGRADVALS(NOISE, MU, VARSIGMA, Y)
%	computes the gradient of the compound noise with respect to the
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
%	% SEEALSO CMPNDNOISEPARAMINIT, CMPNDNOISEGRADIENTPARAM, NOISEGRADVALS, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



startVal = 1;
endVal = 0;
dlnZ_dmu = zeros(size(mu));
dlnZ_dvs = zeros(size(varsigma));
for i = 1:length(noise.comp)
  fhandle = str2func([noise.comp{i}.type 'NoiseGradVals']);
  [dlnZ_dmu(:, i), dlnZ_dvs(:, i)]  = fhandle(noise.comp{i}, ...
                                              mu(:, i), ...
                                              varsigma(:, i), ...
                                              y(:, i));
end
