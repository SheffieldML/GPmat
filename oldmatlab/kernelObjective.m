function f = kernelObjective(parameters, x, h, A, type, prior)

if nargin < 6
  prior = 1;
end
K = kernel(x, parameters, type);
[V, D] = eig(K);
invK = V*diag(1./diag(D))*V';
f = -.5*sum(log(diag(D))) ...
    + .5*trace(eye(length(h))-invK*(h*h'+A));
if prior
  f = f - sum(log(parameters.*parameters));
end
f = -f;
