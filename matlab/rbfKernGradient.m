function g = rbfKernGradient(kern, x, covGrad)

% RBFKERNGRADIENT Gradient of rbf kernel's parameters.

% IVM

[k, dist2xx] = rbfKernCompute(x, kern);
g(1) = - .5*sum(sum(covGrad.*k.*dist2xx))*kern.inverseWidth;
g(2) =  sum(sum(covGrad.*k));
