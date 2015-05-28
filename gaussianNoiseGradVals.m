function [dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y)

% GAUSSIANNOISEGRADVALS Gradient of GAUSSIAN noise log Z with respect to input mean and variance.
% FORMAT
% DESC computes the gradient of the Gaussian
% noise with respect to the input mean and the input variance.
% ARG noise : noise structure for which gradients are being
% computed.
% ARG mu : mean input locations with respect to which gradients are
% being computed.
% ARG varSigma : variance input locations with respect to which
% gradients are being computed.
% ARG y : noise model output observed values associated with the given points.
% RETURN dlnZ_dmu : the gradient of log Z with respect to the input mean.
% RETURN dlnZ_dvs : the gradient of log Z with respect to the input variance.
%
% SEEALSO gaussianNoiseParamInit, gaussianNoiseGradientParam, noiseGradVals, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(y, 2);
nu = 1./(noise.sigma2+varsigma);
dlnZ_dmu = zeros(size(nu));
for i = 1:D
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i) - noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);
