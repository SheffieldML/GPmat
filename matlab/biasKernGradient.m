function g = biasKernGradient(kern, x, covGrad)

% BIASKERNGRADIENT Gradient of bias kernel's parameters.

% KERN

% KERN

g = sum(sum(covGrad));