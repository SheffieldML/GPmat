function L = ivmApproxLogLikelihood(model);

% IVMAPPROXLOGLIKELIHOOD Return the approximate log-likelihood for the IVM.

% IVM

x = model.X(model.I, :);
m = model.m(model.I, :);
K = kernCompute(model.kern, x);
L = 0;

if model.noise.spherical
  % there is only one value for all beta
  [invK, UC] = pdinv(K+diag(1./model.beta(model.I, 1)));
  logDetTerm = logdet(K, UC);
end
  
for i = 1:size(m, 2)
  if ~model.noise.spherical
    [invK, UC] = pdinv(K+diag(1./model.beta(model.I, i)));
    logDetTerm = logdet(K, UC);
  end
  L = L -.5*logDetTerm- .5*m(:, i)'*invK*m(:, i);
end
