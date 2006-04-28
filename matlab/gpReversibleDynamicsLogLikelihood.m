function ll = gpReversibleDynamicsLogLikelihood(model)

% GPREVERSIBLEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.
%
% ll = gpReversibleDynamicsLogLikelihood(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpReversibleDynamicsLogLikelihood.m version 1.1



ll = gpLogLikelihood(model);

% Use prior on first two points only for reversible dynamics model.
if isfield(model, 'prior') &  ~isempty(model.prior)
  ll = ll + priorLogProb(model.prior, model.X(1, :));
  ll = ll + priorLogProb(model.prior, model.X(2, :));
end  
