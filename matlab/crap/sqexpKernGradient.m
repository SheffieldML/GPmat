function g = sqexpKernGradient(kern, x, covGrad)

% SQEXPKERNGRADIENT Gradient of squared exponential kernel's parameters.

% IVM


[k, rbfPart, dist2xx] = sqexpKernCompute(x, kern);
g(1) = - .5*sum(sum(covGrad.*rbfPart.*dist2xx))*kern.inverseWidth;
g(2) =  sum(sum(covGrad.*rbfPart));
g(3) =  sum(sum(covGrad.*eye(size(x, 1))))*kern.whiteVariance;
g(4) =  sum(sum(covGrad))*kern.biasVariance;
