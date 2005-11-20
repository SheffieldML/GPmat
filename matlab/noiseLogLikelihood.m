function L = noiseLogLikelihood(noise, mu, varsigma, y);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.

% NOISE

fhandle = str2func([noise.type 'LogLikelihood']);
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
