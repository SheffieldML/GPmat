function g = rbfardKernGradient(kern, x, covGrad)

% RBFARDKERNGRADIENT Gradient of radial basis function ARD kernel's parameters.

% KERN


g = zeros(1, size(x, 2)+2);
[k, dist2xx] = rbfardKernCompute(kern, x);
covGradK = covGrad.*k;
g(1) = - .5*sum(sum(covGradK.*dist2xx));
g(2) =  sum(sum(covGrad.*k))/kern.variance;

for i = 1:size(x, 2)
  g(2 + i)  =  -(sum(covGradK*(x(:, i).*x(:, i))) ...
                 -x(:, i)'*covGradK*x(:, i))*kern.inverseWidth;
end
