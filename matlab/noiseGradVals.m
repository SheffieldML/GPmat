function [dlnZ_dmu, dlnZ_dvs] = noiseGradVals(noise, mu, varsigma, y)

% NOISEGRADVALS Gradient of noise model wrt mu and varsigma.
% FORMAT
% DESC computes the gradient of the given
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
% SEEALSO noiseCreate, noiseParamInit, noiseGradientParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

fhandle = str2func([noise.type 'NoiseGradVals']);
[dlnZ_dmu, dlnZ_dvs] = fhandle(noise, mu, ...
                               varsigma, y);
%/~
if any(isnan(dlnZ_dmu))
  warning('Gradient of mu is NaN')
end
if any(isnan(dlnZ_dvs))
  warning('Gradient of varsigma is NaN')
end
%~/