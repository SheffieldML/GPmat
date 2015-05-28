function model = mogEstep(model)

% MOGESTEP Do an E-step on an MOG model.
% FORMAT
% DESC carries out an expectation step on a mixtures of
% Gaussians model. 
% ARG model : the model which is to be updated.
% RETURN model : the model with updated posteriors.
%
% SEEALSO : mogCreate, mogUpdateMean, mogUpdateCovariance
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

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
  model.posterior(:, i) = log(model.prior(i))+logDetTerm-0.5* ...
      centredY;
end
maxL = max(model.posterior, [], 2);
model.posterior = model.posterior - repmat(maxL, 1, model.m);
model.lnposterior = model.posterior;
model.posterior = exp(model.posterior);
sumPosterior = repmat(sum(model.posterior, 2), 1, model.m);
model.posterior = model.posterior./sumPosterior;
model.lnposterior = model.lnposterior - log(sumPosterior);
%/~
if any(any(isnan(model.posterior)))
  error('Posterior is nan');
end
%~/
