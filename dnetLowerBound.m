function lll = dnetLowerBound(model)

% DNETLOWERBOUND Computes lower bound on log likelihood for an DNET model.
% FORMAT
% DESC computes the variational lower bound on the log likelihood
% of a mixtures of probabilistic PCA model.
% ARG model : the model for which log likelihood is to be computed.
% RETURN lll : the lower bound on the log likelihood computed for the model.
% 
% SEEALSO : dnetCreate, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS
  
lll = 0.5*model.N*model.d*log(model.beta/(2*pi)) ...
      - model.N*log(model.M);

lll = lll - 0.5*model.alpha*sum(sum(model.A.*model.A)) ...
      - 0.5*model.alpha*sum(model.b.*model.b);

lll = lll - sum(sum(xlogy(model.w)));


% Get projections of latent samples.
Ypred = dnetOut(model);
if model.N > model.M
  for i = 1:model.M
    diffY = model.y - repmat(Ypred(i, :), model.N, 1);
    diffY = diffY.*diffY.*repmat(model.w(:,i), 1, model.d);
    lll = lll - 0.5*model.beta*sum(sum(diffY));
  end
else
  for i = 1:model.N
    diffY = repmat(model.y(i, :), model.M, 1) - Ypred;
    diffY =diffY.*diffY.*repmat(model.w(i,:)', 1, model.d);
    lll = lll - 0.5*model.beta*sum(sum(diffY));
  end
end

