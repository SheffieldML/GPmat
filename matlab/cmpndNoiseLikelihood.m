function L = cmpndNoiseLikelihood(noise, mu, varsigma, y)


% CMPNDNOISELIKELIHOOD Likelihood of the data under the CMPND noise model.
% FORMAT
% DESC returns the likelihood of a data set under the  compound noise model.
% ARG noise : the noise structure for which the likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : cmpndNoiseParamInit, cmpndNoiseLogLikelihood, noiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


L = zeros(size(mu));
for i = 1:length(noise.comp)
  L(:, i) = noiseLikelihood(noise.comp{i},...
                mu(:, i), ...
                varsigma(:, i), ...
                y(:, i));
end