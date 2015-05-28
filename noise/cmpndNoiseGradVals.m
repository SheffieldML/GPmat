function [dlnZ_dmu, dlnZ_dvs] = cmpndNoiseGradVals(noise, mu, varsigma, y)


% CMPNDNOISEGRADVALS Gradient of CMPND noise log Z with respect to input mean and variance.
% FORMAT
% DESC computes the gradient of the compound
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
% SEEALSO cmpndNoiseParamInit, cmpndNoiseGradientParam, noiseGradVals, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


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
