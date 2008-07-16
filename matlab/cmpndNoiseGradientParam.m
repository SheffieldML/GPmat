function g = cmpndNoiseGradientParam(noise, mu, varsigma, y)


% CMPNDNOISEGRADIENTPARAM Gradient of CMPND noise's parameters.
% FORMAT
% DESC computes the gradient of the log Z of the compound noise model with respect to the of functions with respect to the
% compound
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
% SEEALSO cmpndNoiseParamInit, cmpndnoiseGradVals, noiseGradientParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


g = zeros(1, noise.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  g(1, startVal:endVal)  = noiseGradientParam(noise.comp{i}, ...
					      mu(:, i), ...
					      varsigma(:, i), ...
					      y(:, i));
  startVal = endVal + 1;
end
g = g*noise.paramGroups;