function g = rbfKernGradient(kern, x, covGrad)

% RBFKERNGRADIENT Gradient of rbf kernel's parameters.

% IVM

[k, dist2xx] = rbfKernCompute(kern, x);
g(1) = - .5*sum(sum(covGrad.*k.*dist2xx));
g(2) =  sum(sum(covGrad.*k))/kern.variance;

%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/