function g = rbfKernGradient(kern, x, covGrad)

% RBFKERNGRADIENT Gradient of rbf kernel's parameters.

% IVM

[k, dist2xx] = rbfKernCompute(kern, x);
g(1) = - .5*sum(sum(covGrad.*k.*dist2xx));
g(2) =  sum(sum(covGrad.*k))/kern.variance;

if kern.linearBound
  g(1) = g(1)*gradFactLinearBound(kern.inverseWidth);
  g(2) = g(2)*gradFactLinearBound(kern.variance);
else
  g(1) = g(1)*kern.inverseWidth;
  g(2) = g(2)*kern.variance;
end

%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/