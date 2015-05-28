function ll = dnetLogLikelihood(model)

% DNETLOGLIKELIHOOD Density network log likelihood.
% FORMAT
% DESC computes the log likelihood of a density network
% model. 
% ARG model : the model structure for computing the log likelihood.
% RETURN ll : the model log likelihood.
%
% SEEALSO : modelLogLikeihood
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

  
ll = 0.5*model.d*model.N*log(model.beta/(2*pi)) ...
     - model.N*log(model.M); 

ll = ll - 0.5*model.alpha*sum(sum(model.A.*model.A)) ...
     - 0.5*model.alpha*sum(model.b.*model.b);

  
llPointComp = zeros(model.N, model.M);
% Get projections of latent samples.
Ypred = dnetOut(model);

if model.N > model.M
  for i = 1:model.M
    diffY = model.y - repmat(Ypred(i, :), model.N, 1);
    diffY = diffY.*diffY;
    llPointComp(:, i) = - 0.5*model.beta*sum(diffY, 2);
  end
else
  for i = 1:model.N
    diffY = repmat(model.y(i, :), model.M, 1) - Ypred;
    diffY = diffY.*diffY;
    llPointComp(i, :) = - 0.5*model.beta*sum(diffY, 2)';
  end
end

maxllPointComp = max(llPointComp, [], 2);
llPointComp = exp(llPointComp - repmat(maxllPointComp, 1, model.M));
ll = ll + sum(log(sum(llPointComp, 2)) + maxllPointComp);
