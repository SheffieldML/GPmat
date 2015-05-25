function g = fgplvmSequenceLogLikeGradient(model, X, Y, varargin)

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
% ARG P1, P2, P3 ... : optional additional arguments to be passed
% to the dynamics' model sequence log likelihood gradient.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent sequence.
%
% SEEALSO : fgplvmSequenceLogLikelihood, fgplvmOptimiseSequence, fgplvmSequenceLogLikeGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Carl Henrik Ek, 2008

% FGPLVM
if(nargin<3)
  error('This function requires at least two arguments');
end

g = gpSequenceLogLikeGradient(model, X, Y);

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  % A dynamics model is being used.
  feval = str2func([model.dynamics.type ...
		    'SequenceLogLikeGradient']);
  if(isfield(model.dynamics,'indexIn')&&~isempty(model.dynamics.indexIn))
    dim = model.dynamics.indexIn;
  else
    dim = 1:1:size(g,2);
  end
  % 'balancing' of the dynamics alla Urtasun.
  if isfield(model, 'balancing') & ~isempty(model.balancing)
    g(:,dim) = g(:,dim) + model.balancing*feval(model.dynamics, X);
  else
    g = g + feval(model.dynamics, X);
  end
  g = g + feval(model.dynamics, X, varargin{:});
elseif isfield(model, 'prior') &  ~isempty(model.prior)
  for i = 1:size(X, 1)
    g(i, :) = g(i, :) + priorGradient(model.prior, X(i, :));
  end
end
  

  
