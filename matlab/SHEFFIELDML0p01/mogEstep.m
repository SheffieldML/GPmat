function model = mogEstep(model)

% MOGESTEP Do an E-step on an MOG model.
%
%	Description:
%
%	MODEL = MOGESTEP(MODEL) carries out an expectation step on a
%	mixtures of Gaussians model.
%	 Returns:
%	  MODEL - the model with updated posteriors.
%	 Arguments:
%	  MODEL - the model which is to be updated.
%	
%
%	See also
%	MOGCREATE, MOGUPDATEMEAN, MOGUPDATECOVARIANCE


%	Copyright (c) 2006 Neil D. Lawrence


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
