function g = kernelGradient(params, model, prior)

% KERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.

% IVM
%/~
if any(isnan(params))
  warning('Parameter is NaN')
end
%~/
if nargin < 3
  prior = 1;
end

x = model.X(model.I, :);
m = model.m(model.I, :);
model.kern = kernExpandParam(model.kern, params);
K = kernCompute(model.kern, x);
g = zeros(size(params));
for j = 1:size(m, 2)
  invK = pdinv(K+diag(1./model.beta(model.I, j)));
  covGrad = covarianceGradient(invK, m(:, j));
  g = g + kernGradient(model.kern, x, covGrad);
end  
if prior
  g = g - 1;
end
g = -g;

