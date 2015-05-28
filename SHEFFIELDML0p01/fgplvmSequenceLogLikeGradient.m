function g = fgplvmSequenceLogLikeGradient(model, X, Y, varargin)

% FGPLVMSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM.
%
%	Description:
%
%	G = FGPLVMSEQUENCELOGLIKEGRADIENT(MODEL, X, Y, ...) returns the
%	gradient of the log likelihood with respect to the latent position,
%	where the log likelihood is conditioned on the training set.
%	 Returns:
%	  G - the gradient of the log likelihood, conditioned on the
%	   training data, with respect to the latent sequence.
%	 Arguments:
%	  MODEL - the model for which the gradient computation is being
%	   done.
%	  X - the latent sequence where the gradient is being computed.
%	  Y - the sequence in data space for which the computation is being
%	   done.
%	  ... - optional additional arguments to be passed to the dynamics'
%	   model sequence log likelihood gradient.
%	
%	
%
%	See also
%	FGPLVMSEQUENCELOGLIKELIHOOD, FGPLVMOPTIMISESEQUENCE, FGPLVMSEQUENCELOGLIKEGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence


%	With modifications by Carl Henrik Ek 2008

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
  

  
