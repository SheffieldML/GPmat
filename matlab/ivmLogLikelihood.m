function L = ivmLogLikelihood(model, x, y);

% IVMLOGLIKELIHOOD Return the log-likelihood for the IVM.

% IVM

if nargin < 3
  % This implies evaluate for the training data.
  mu = model.mu;
  varsigma = model.varSigma;
  y = model.y;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

L = noiseLogLikelihood(model.noise, mu, varsigma, y);

% check if there is a prior over kernel parameters
if isfield(model.kern, 'priors')
  params = feval([kern.type 'KernExpandParams'], model.kern);
  for i = 1:length(model.kern.priors)
    index = model.kern.priors(i).index;
    L = L + priorLogProb(model.kern.priors(i), params(index));
  end
end
