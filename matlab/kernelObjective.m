function f = kernelObjective(params, model, prior)

% KERNELOBJECTIVE Likelihood approximation.

% KERN

%/~
if any(isnan(params))
  warning('Parameter is NaN')
end
%~/

if nargin < 3
  prior = 0;
end
x = model.X(model.I, :);
m = model.m(model.I, :);
model.kern = kernExpandParam(model.kern, params);
K = kernCompute(model.kern, x);
f = 0;

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
  f = f -.5*logDetTerm- .5*m(:, i)'*invK*m(:, i);
end
f = f + kernPriorLogProb(model.kern);


f = -f;
