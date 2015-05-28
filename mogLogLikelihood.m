function ll = mogLogLikelihood(model)
  
% MOGLOGLIKELIHOOD Mixture of Gaussian's log likelihood.
% FORMAT
% DESC computes the variational lower bound on the log likelihood
% of a mixtures of probabilistic PCA model, it wraps the mogLowerBound command.
% ARG model : the model for which log likelihood is to be computed.
% RETURN lll : the lower bound on the log likelihood computed for the model.
% 
% SEEALSO : mogCreate, mogLowerBound
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008
% MLTOOLS
  
model = mogEstep(model);
ll = mogLowerBound(model);

