function f = kernelObjective(params, model, prior)

% KERNELOBJECTIVE Likelihood approximation.

% IVM

if nargin < 3
  prior = 1;
end
x = model.X(model.I, :);
m = model.m(model.I, :);
model.kern = kernExpandParam(params, model.kern);
K = kernCompute(x, model.kern);
f = 0;
for i = 1:size(m, 2)
  [invK, UC] = pdinv(K+diag(1./model.beta(model.I, i)));
  f = f -.5*logdet(K, UC) - .5*m(:, i)'*invK*m(:, i);
end
if prior
  f = f - sum(params);
end
f = -f;
