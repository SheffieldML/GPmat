function [dlnZ_dmu, dlnZ_dvs] = probitNoiseGradVals(noise, mu, varsigma, y)


% PROBITNOISEGRADVALS Gradient of PROBIT noise log Z with respect to input mean and variance.
% FORMAT
% DESC computes the gradient of the probit based classification
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
% SEEALSO probitNoiseParamInit, probitNoiseGradientParam, noiseGradVals, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(mu, 2);
c = y./sqrt(noise.sigma2+varsigma);
u = zeros(size(c));
for i = 1:D
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
end
dlnZ_dmu = c.*gradLogCumGaussian(u);
dlnZ_dvs = -.5*c.*u.*dlnZ_dmu;
