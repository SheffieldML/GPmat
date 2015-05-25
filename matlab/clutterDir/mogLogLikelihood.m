function ll = mogLogLikelihood(model)

% MOGLOGLIKELIHOOD Computes log likelihood for an MOG model.
% FORMAT
% DESC computes the log likelihood of a mixtures of Gaussians
% model.
% ARG model : the model for which log likelihood is to be computed.
% RETURN ll : the log likelihood computed for the model.
% 
% SEEALSO : mogCreate, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006

for i = 1:model.m
  centredY = model.Y - repmat(model.mean(i, :), model.N, 1);
  switch model.covtype
    case 'ppca'
     centredY = centredY/model.U{i};
     logDetTerm = -0.5*logdet([], model.U{i});
   case 'spherical'
    centredY = centredY/sqrt(model.sigma2(i));
    logDetTerm = -0.5*model.d*log(model.sigma2(i));
  end
  centredY = sum(centredY.*centredY, 2);
  L(:, i) = logDetTerm-0.5* ...
            centredY-model.d/2*log(2*pi) + log(model.prior(i));
end
maxL = max(L, [], 2);
L = L - repmat(maxL, 1, model.m);
L = log(sum(exp(L), 2))+maxL;
ll = sum(L);