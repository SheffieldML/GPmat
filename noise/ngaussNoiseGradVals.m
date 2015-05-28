function [dlnZ_dmu, dlnZ_dvs] = ngaussNoiseGradVals(noise, mu, varsigma, y)


% NGAUSSNOISEGRADVALS Gradient of NGAUSS noise log Z with respect to input mean and variance.
% FORMAT
% DESC computes the gradient of the noiseless Gaussian
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
% SEEALSO ngaussNoiseParamInit, ngaussNoiseGradientParam, noiseGradVals, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


[dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y);