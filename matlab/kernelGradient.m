function g = kernelGradient(params, model, prior)

% KERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.

% IVM

if nargin < 3
  prior = 1;
end

x = model.X(model.I, :);
m = model.m(model.I, :);
params = log(thetaConstrain(exp(params)));
model.kern = kernExpandParam(params, model.kern);
K = kernCompute(x, model.kern);
g = zeros(size(params));
for j = 1:size(m, 2)
  invK = pdinv(K+diag(1./model.beta(model.I, j)));
  covGrad = covarianceGradient(invK, m(:, j));
  g = g + feval([model.kern.type 'KernGradient'], model.kern, x, covGrad);
end  
if prior
  g = g - 1;
end
g = -g;

