function L = cmpndNoiseLogLikelihood(noise, mu, varsigma, y)


% CMPNDNOISELOGLIKELIHOOD Log likelihood of the data under the CMPND noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the  compound noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : cmpndNoiseParamInit, cmpndNoiseLikelihood, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


L = 0;
for i = 1:length(noise.comp)
  L = L + noiseLogLikelihood(noise.comp{i}, ...
			     mu(:, i), ...
			     varsigma(:, i), ...
			     y(:, i));
end