function g = ivmApproxLogLikeKernGrad(model)

% IVMAPPROXLOGLIKEKERNGRAD Gradient of the approximate likelihood wrt kernel parameters.

% IVM

x = model.X(model.I, :);
m = model.m(model.I, :);
K = kernCompute(model.kern, x);
g = zeros(1, model.kern.nParams);

if model.noise.spherical
  % there is only one value for all beta
  invK = pdinv(K+diag(1./model.beta(model.I, 1)));
end

for j = 1:size(m, 2)
  if ~model.noise.spherical
    invK = pdinv(K+diag(1./model.beta(model.I, j)));
  end
  covGrad = feval([model.type 'CovarianceGradient'], invK, m(:, j));
  g = g + kernGradient(model.kern, x, covGrad);
end  
