function L = noiseLogLikelihood(noise, mu, varsigma, y);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.

% IVM

L = feval([noise.type 'LogLikelihood'], noise, mu, varsigma, y);


% check if there is a prior over parameters
if isfield(noise, 'priors')
  params = feval([noise.type 'NoiseExpandParams'], noise);
  for i = 1:length(noise.priors)
    index = noise.priors(i).index;
    L = L + priorLogProb(noise.priors(i), params(index));
  end
end
