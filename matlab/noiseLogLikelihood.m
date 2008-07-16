function L = noiseLogLikelihood(noise, mu, varsigma, y);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the given noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : noiseParamInit, noiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

fhandle = str2func([noise.type 'NoiseLogLikelihood']);
L = fhandle(noise, mu, varsigma, y);


% check if there is a prior over parameters
if isfield(noise, 'priors')
  fhandle = str2func([noise.type 'NoiseExpandParams']);
  params = fhandle(noise);
  for i = 1:length(noise.priors)
    index = noise.priors(i).index;
    L = L + priorLogProb(noise.priors(i), params(index));
  end
end
