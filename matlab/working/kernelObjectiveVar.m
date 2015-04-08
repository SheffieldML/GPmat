function f = kernelObjectiveVar(parameters, x, h, A, type, prior)

% KERNELOBJECTIVEVAR Variational bound for kernel parameters.

if nargin < 6
  prior = 1;
end
parameters=log(thetaConstrain(exp(parameters)));

K = kernel(x, parameters, type);
%[V, D] = eig(K);
%invK = V*diag(1./diag(D))*V';

invK = pdinv(K);
f = -.5*logdet(K) ...
    + .5*trace(eye(size(h, 1))-invK*(h*h'+A));
if prior
  f = f - sum(parameters);
end
f = -f;
