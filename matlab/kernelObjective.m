function f = kernelObjective(params, model, prior)

% KERNELOBJECTIVE Likelihood approximation.

% IVM
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

if strcmp(model.noise.type, 'gaussian')
  [invK, UC] = pdinv(K+diag(1./model.beta(model.I, 1)));
  logDetTerm = logDet(K, UC);
end
  
for i = 1:size(m, 2)
  if ~strcmp(model.noise.type, 'gaussian')
    [invK, UC] = pdinv(K+diag(1./model.beta(model.I, i)));
    logDetTerm = logDet(K, UC);
  end
    f = f -.5*logDetTerm- .5*m(:, i)'*invK*m(:, i);
end
if prior
  f = f - sum(params);
end
f = -f;
