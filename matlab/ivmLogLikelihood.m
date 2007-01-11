function L = ivmLogLikelihood(model, x, y);

% IVMLOGLIKELIHOOD Return the log-likelihood for the IVM.
% FORMAT
% DESC computes the log likelihood of the given data points under
% the IVM and the given noise model.
% ARG model : the model for which the log likelihood is computed.
% ARG x : the input locations.
% ARG y : the target values.
%
% SEEALSO : ivmNegLogLikelihood, ivmPosteriorMeanVar,
% noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence

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
  fhandle = str2func([model.kern.type 'KernExpandParams']);
  params = fhandle(model.kern);
  for i = 1:length(model.kern.priors)
    index = model.kern.priors(i).index;
    L = L + priorLogProb(model.kern.priors(i), params(index));
  end
end
