function lll = dnetLowerBound(model)

% DNETLOWERBOUND Computes lower bound on log likelihood for an DNET model.
%
%	Description:
%
%	LLL = DNETLOWERBOUND(MODEL) computes the variational lower bound on
%	the log likelihood of a mixtures of probabilistic PCA model.
%	 Returns:
%	  LLL - the lower bound on the log likelihood computed for the
%	   model.
%	 Arguments:
%	  MODEL - the model for which log likelihood is to be computed.
%	
%
%	See also
%	DNETCREATE, MODELLOGLIKELIHOOD


%	Copyright (c) 2008 Neil D. Lawrence

  
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

