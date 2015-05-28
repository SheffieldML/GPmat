function [dlnZ_dmu, dlnZ_dvs] = mgaussianNoiseGradVals(noise, mu, varsigma, y)


% MGAUSSIANNOISEGRADVALS Gradient of MGAUSSIAN noise log Z with respect to input mean and variance.
% FORMAT
% DESC computes the gradient of the multiple output Gaussian
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
% SEEALSO mgaussianNoiseParamInit, mgaussianNoiseGradientParam, noiseGradVals, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


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
%/~
if any(isnan(dlnZ_dmu))
  warning('dlnZ_dmu is NaN')
end
if any(isnan(dlnZ_dvs))
  warning('dlnZ_dvs is NaN')
end
%~/