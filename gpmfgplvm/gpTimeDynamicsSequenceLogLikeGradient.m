function gX = gpTimeDynamicsSequenceLogLikeGradient(model, latentValsRaw, t)

% GPTIMEDYNAMICSSEQUENCELOGLIKEGRADIENT Log-likelihood gradient for of a sequence of the GP-LVM time dynamics.
% FORMAT
% DESC returns the gradient of the log likelihood with respect to
% the latent position, where the log likelihood is conditioned on
% the training set latent points. 
% ARG model : the model for which the gradient computation is being
% done.
% ARG X : the latent sequence where the gradient is being computed.
% ARG t : the times in the sequence for which the log
% likelihood is to be computed.
% RETURN g : the gradient of the log likelihood, conditioned on the
% training data, with respect to the latent sequence.
%
% SEEALSO : gpDynamicsSequenceLogLikelihood, fgplvmOptimiseSequence, fgplvmSequenceLogLikeGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Carl Henrik Ek, 2008
%

% FGPLVM

if nargin < 3
  t = [2:size(latentValsRaw, 1)]';
  warning('No time inputs provided for gpTimeDynamicsSequenceLogLikelihood')
end
if(isfield(model,'indexIn')&&~isempty(model.indexIn))
  latentValsRaw = latentValsRaw(:,model.indexIn);
end

ind_in = 1:size(latentValsRaw, 1)-1;
ind_out = 2:size(latentValsRaw, 1);

latentVals = latentValsRaw(ind_in, :);
if model.diff
  Y = latentValsRaw(ind_out, :) - latentVals;
else
  Y = latentValsRaw(ind_out, :);
end

gX = zeros(size(latentValsRaw, 1), size(latentValsRaw, 2));

[mu, covar, factors] = gpPosteriorMeanCovar(model, t);
invCovar = pdinv(covar);
Ydiff= Y - mu;
% Deal with the fact that X appears in the *target* for the dynamics.
for i =1:model.d
  gX(ind_out, i) = gX(ind_out, i) - 1/(factors(i))*invCovar*Ydiff(:, i);
  if model.diff
    gX(ind_in, i) = gX(ind_in, i) + 1/(factors(i))*invCovar*Ydiff(:, i);
  end
end
 


gX(ind_in,:) = gX(ind_in,:);

%/~
% Maybe need to consider using a prior here.
% if isfield(model, 'prior') &  ~isempty(model.prior)
%   gX(1,:) = gX(1,:) + priorGradient(model.prior, model.X(1,:));
% end
%~/
