function lll = mogLowerBound(model)

% MOGLOWERBOUND Computes lower bound on log likelihood for an MOG model.
% FORMAT
% DESC computes the variational lower bound on the log likelihood
% of a mixtures of probabilistic PCA model.
% ARG model : the model for which log likelihood is to be computed.
% RETURN lll : the lower bound on the log likelihood computed for the model.
% 
% SEEALSO : mogCreate, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008

% MLTOOLS

lll = -sum(sum(xlogy(model.posterior)));

if model.isInfinite
  % DP add in extra posterior entropy etc.
end

for i = 1:model.m
  centredY = model.Y - repmat(model.mean(i, :), model.N, 1);
  switch model.covtype
   case 'ppca'
    centredY = centredY/model.U{i};
    logDetTerm = logdet([], model.U{i});
   case 'spherical'
    centredY = centredY/sqrt(model.sigma2(i));
    logDetTerm = model.d*log(model.sigma2(i));
  end
  centredY = sum(centredY.*centredY, 2);
  lll = lll - 0.5*sum(model.posterior(:, i).*(logDetTerm+centredY-2*log(model.prior(i)+1e-300)));
end
lll = lll - model.N*model.d/2*log(2*pi);  
