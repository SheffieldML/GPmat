function ll = gpDynamicsLogLikelihood(model)

% GPDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.
%
% ll = gpDynamicsLogLikelihood(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpDynamicsLogLikelihood.m version 1.1



ll = gpLogLikelihood(model);

% Use prior on first point only for dynamics models.
if isfield(model, 'prior') &  ~isempty(model.prior)
  ll = ll + priorLogProb(model.prior, model.X(1, :));
end  
