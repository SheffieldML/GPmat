function g = rbfardKernGradient(kern, x, covGrad)

% RBFARDKERNGRADIENT Gradient of radial basis function ARD kernel's parameters.

% IVM

g = zeros(1, size(x, 2)+2);
[k, dist2xx] = rbfardKernCompute(kern, x);
covGradK = covGrad.*k;
g(1) = - .5*sum(sum(covGradK.*dist2xx));
g(2) =  sum(sum(covGrad.*k))/kern.variance;

for i = 1:size(x, 2)
  g(2 + i)  =  -(sum(covGradK*(x(:, i).*x(:, i))) ...
                 -x(:, i)'*covGradK*x(:, i)) ...
      *kern.inputScales(i)*(1-kern.inputScales(i))*kern.inverseWidth;
end

if kern.linearBound
  g(1) = g(1) * gradFactLinearBound(kern.inverseWidth);
  g(2) = g(2) * gradFactLinearBound(kern.variance);
else
  g(1) = g(1) *kern.inverseWidth;
  g(2) = g(2) *kern.variance;
end