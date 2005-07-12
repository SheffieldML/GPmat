function g = sqexpKernGradient(kern, x, covGrad)

% SQEXPKERNGRADIENT Gradient of squared exponential kernel's parameters.

% KERN

[k, rbfPart, dist2xx] = sqexpKernCompute(kern, x);
g(1) = - .5*sum(sum(covGrad.*rbfPart.*dist2xx));
g(2) =  sum(sum(covGrad.*rbfPart))/kern.rbfVariance;
g(3) =  sum(sum(covGrad));
g(4) =  trace(covGrad);

