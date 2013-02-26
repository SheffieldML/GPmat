function ll = gpDynamicsSequenceLogLikelihood(model, Xraw, varargin)

% GPDYNAMICSSEQUENCELOGLIKELIHOOD Return the log likelihood of a given latent sequence.
%
%	Description:
%
%	LL = GPDYNAMICSSEQUENCELOGLIKELIHOOD(MODEL, X) returns a log
%	likelihood for a latent sequence given a dynamics model and a set of
%	training data.
%	 Returns:
%	  LL - the log probability of the latent sequence given the training
%	   data and the model.
%	 Arguments:
%	  MODEL - the model containing the dynamics model for which the log
%	   likelihood is to be computed.
%	  X - the latent positions of the sequence for which the log
%	   likelihood is to be computed.
%	
%	
%
%	See also
%	FGPLVMSEQUENCELOGLIKELIHOOD, GPDYNAMICSCREATE


%	Copyright (c) 2006 Neil D. Lawrence


%	With modifications by Carl Henrik Ek 2008


if(isfield(model,'indexIn')&&~isempty(model.indexIn)&&length(model.indexIn)~=size(Xraw,2))
  Xraw = Xraw(:,model.indexIn);
end

startVal=1;
X = Xraw(1:end-1, :);

if model.diff
  Y = Xraw(2:end, :) - X;
else
  Y = Xraw(2:end, :);
end

logTwoPi = log(2*pi);
[mu, covarSigma, factors] = gpPosteriorMeanCovar(model, X);
U = jitChol(covarSigma);
Ydiff = (Y-mu);
ll =0;
UinvYdiff = U'\Ydiff;
logDet = logdet([], U);
for i = 1:model.d
  ll = ll + logDet + log(factors(i)) + (UinvYdiff(:, i)'*UinvYdiff(:, ...
                                                    i))/factors(i);
  ll = ll + logTwoPi*size(Y, 1);
end
ll = -0.5*ll;

% Check for prior on first point
if isfield(model, 'prior') & ~isempty(model.prior)
  ll = ll + priorLogProb(model.prior, X(1, :));
end