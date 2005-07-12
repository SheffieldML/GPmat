function g = biasKernGradient(kern, x, covGrad)

% BIASKERNGRADIENT Gradient of bias kernel's parameters.

% KERN

g = sum(sum(covGrad));
