function ll = fgplvmSequenceLogLikelihood(model, X, Y)

% FGPLVMSEQUENCELOGLIKELIHOOD Log-likelihood of a sequence for the GP-LVM.
% FORMAT
% DESC returns the log probability of a given latent sequence and an
% associated observed data sequence under the posterior prediction
% induced by the training data for a given GP-LVM model.
% ARG model : the model for which the sequence prediction will be
% made.
% ARG X : the latent sequence for which the posterior distribution
% will be evaluated.
% ARG Y : the observed data sequence for which the posterior
% distribution will be evaluated.
%
% SEEALSO : fgplvmCreate, fgplvmOptimiseSequence, fgplvmSequenceObjective
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% FGPLVM

logTwoPi = log(2*pi);
[mu, covarSigma, factors] = gpPosteriorMeanCovar(model, X);
missing = true;
if ~any(isnan(Y))
  U = jitChol(covarSigma);
  missing = false;
end
Ydiff = (Y-mu);
ll =0;
for i = 1:model.d
  if missing
    ind = find(~isnan(Ydiff(:, i)));
    if length(ind) ~= 0    
      U = jitChol(covarSigma(ind, ind));
    end
  else
    ind = [1:size(Ydiff, 1)]';
  end
  if length(ind) ~= 0    
    UinvYdiff = U'\Ydiff(ind, i);
    logDet = logdet([], U);
    ll = ll + logDet + log(factors(i)) + (UinvYdiff'*UinvYdiff)/factors(i);
    ll = ll + logTwoPi*length(ind);
  end
end
ll = -0.5*ll;


if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  % A dynamics model is being used.
  feval = str2func([model.dynamics.type 'SequenceLogLikelihood']);
  ll = ll + feval(model.dynamics, X);
elseif isfield(model, 'prior') &  ~isempty(model.prior)
  for i = 1:size(X, 1)
    ll = ll + priorLogProb(model.prior, X(i, :));
  end
end
