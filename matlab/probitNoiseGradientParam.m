function g = probitNoiseGradientParam(noise, mu, varsigma, y)


% PROBITNOISEGRADIENTPARAM Gradient of PROBIT noise's parameters.
% FORMAT
% DESC computes the gradient of the log Z of the probit based classification noise model with respect to the of functions with respect to the
% probit based classification
% noise's parameters. 
% ARG noise : the noise structure for which the gradients are being
% computed.
% ARG mu : the input means for which the gradients are being computed.
% ARG varSigma : the input variances for which the gradients are being computed.
% ARG y : the target values for the noise model.
% RETURN g : gradients of the log Z with respect to
% the noise parameters. The ordering of the vector should match
% that provided by the function noiseExtractParam.
%
%
% SEEALSO probitNoiseParamInit, probitnoiseGradVals, noiseGradientParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


c = y./sqrt(noise.sigma2 + varsigma);
for i = 1:size(mu, 2)
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
end
g = sum(c.*gradLogCumGaussian(u), 1);