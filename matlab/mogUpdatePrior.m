function model = mogUpdatePrior(model)

% MOGUPDATEPRIORML Update the priors of an MOG model.
% FORMAT
% DESC updates the prior probabilities of a mixtures of
% Gaussians model. 
% ARG model : the model which is to be updated.
% RETURN model : the model with updated priors.
%
% SEEALSO : mogCreate, mogUpdateMean, mogUpdateCovariance, mogEstep
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008

% MLTOOLS

if model.isInfinite
  % First compute expectations of v.
  sumS = sum(model.posterior);
  a0bar = model.a0 + sumS; % Posterior value for a_0.
  a1bar = model.a1 + cumsum(sumS); % Posterior value for a_1.
  model.v = a0bar./(a0bar+a1bar);
  tmp = cumprod(1-model.v);
  model.prior = model.v;
  model.prior(2:end) = model.prior(2:end).*tmp(1:end-1);
else
  model.prior = mean(model.posterior);
  model.prior(find(model.prior==0))=1e-100;
end
