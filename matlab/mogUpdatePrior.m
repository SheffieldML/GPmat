function model = mogUpdatePrior(model)

% MOGUPDATEPRIOR Update the priors of an MOG model.
% FORMAT
% DESC updates the prior probabilities of a mixtures of
% Gaussians model. 
% ARG model : the model which is to be updated.
% RETURN model : the model with updated priors.
%
% SEEALSO : mogCreate, mogUpdateMean, mogUpdateCovariance, mogEstep
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

model.prior = mean(model.posterior);

