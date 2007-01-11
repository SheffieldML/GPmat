function L = ivmLikelihoods(model, x, y);

% IVMLOGLIKELIHOODS Return the likelihood for each point for the IVM.
% FORMAT
% DESC computes the likelihood of each given point under the given
% IVM model structure.
% ARG model : the IVM for which the likelihoods are to be computed.
% ARG x : the input points where the likelihoods are to be
% computed.
% ARG y : the target points where the likelihoods are to be
% computed.
%
% SEEALSO : noiseLikelihood, ivmPosteriorMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2005

if nargin < 3
  % This implies evaluate for the traing data.
  mu = model.mu;
  varsigma = model.varSigma;
  y = model.y;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

L = noiseLikelihood(model.noise, mu, varsigma, y);
