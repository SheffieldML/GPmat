function g = fgplvmSequenceLogLikeGradient(model, X, Y)

% FGPLVMSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM.
% FORMAT
% DESC returns the gradient of the log likelihood with respect to
% the latent position, where the log likelihood is conditioned on
% the training set. 
% ARG model : the model for which the gradient computation is being
% done.
% ARG X : the latent sequence where the gradient is being computed.
% ARG Y : the sequence in data space for which the computation is
% being done.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent sequence.
%
% SEEALSO : fgplvmSequenceLogLikelihood, fgplvmOptimiseSequence, fgplvmSequenceLogLikeGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

g = gpSequenceLogLikeGradient(model, X, Y);

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  % A dynamics model is being used.
  feval = str2func([model.dynamics.type 'SequenceLogLikeGradient']);
  g = g + feval(model.dynamics, X);
elseif isfield(model, 'prior') &  ~isempty(model.prior)
  for i = 1:size(X, 1)
    g(i, :) = g(i, :) + priorGradient(model.prior, X(i, :));
  end
end
  

  
