function g = pmvuLogLikeGradients(model)

% PMVULOGLIKEGRADIENTS Gradient of PMVU model log likelihood with respect to parameters.
% FORMAT
% DESC computes the gradient of the probabilistic maximum variance unfolding
% model's log likelihood with respect to the parameters.
% ARG model : model structure for which gradients are being
% computed.
% RETURN g : the returned gradients. 
%
% SEEALSO pmvuCreate, pmvuLogLikelihood, modelLogLikeGradients 
%
% COPYRIGHT : Neil D. Lawrence 2009

% MLTOOLS

expD2 = zeros(model.N, model.k);
for i = 1:model.N
  for j = 1:model.k
    ind = model.indices(i, j);
    D2_e(i, j) = model.K(i, i) - 2*model.K(i, ind) + model.K(ind, ind);
  end
end
gV = 0.5*( model.d*D2_e - model.D2);
fhandle = str2func([model.kappaTransform 'Transform']);
factors = fhandle(model.kappa, 'gradfact');
gV = gV.*factors;
g = gV(:)';
