function g = rbfardKernGradient(kern, x, covGrad)

% RBFARDKERNGRADIENT Gradient of radial basis function ARD kernel's parameters.

% IVM

g = zeros(1, size(x, 2)+2);
[k, dist2xx] = rbfardKernCompute(x, kern);
covGradK = covGrad.*k;
g(1) = - .5*sum(sum(covGradK.*dist2xx))*kern.inverseWidth;
g(2) =  sum(sum(covGrad.*k));

for i = 1:size(x, 2)
  g(2 + i)  =  -(sum(covGradK*(x(:, i).*x(:, i))) ...
                 -x(:, i)'*covGradK*x(:, i)) ...
      *kern.inputScales(i)*(1-kern.inputScales(i))*kern.inverseWidth;
end