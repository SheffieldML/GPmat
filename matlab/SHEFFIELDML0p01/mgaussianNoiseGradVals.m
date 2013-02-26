function [dlnZ_dmu, dlnZ_dvs] = mgaussianNoiseGradVals(noise, mu, varsigma, y)

% MGAUSSIANNOISEGRADVALS Gradient of MGAUSSIAN noise log Z with respect to input mean and variance.
%
%	Description:
%
%	[DLNZ_DMU, DLNZ_DVS] = MGAUSSIANNOISEGRADVALS(NOISE, MU, VARSIGMA,
%	Y) computes the gradient of the multiple output Gaussian noise with
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
%	% SEEALSO MGAUSSIANNOISEPARAMINIT, MGAUSSIANNOISEGRADIENTPARAM, NOISEGRADVALS, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



D = size(y, 2);
nu = zeros(size(y));
dlnZ_dmu = zeros(size(y));
for i = 1:D
  nu(:, i) = 1./(noise.sigma2(i)+varsigma(:, i));
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i) - noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = -.5*nu+.5*dlnZ_dmu.*dlnZ_dmu;

% Remove missing values.
dlnZ_dmu(find(isnan(y))) = 0;
dlnZ_dvs(find(isnan(y))) = 0;
