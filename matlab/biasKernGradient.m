function g = biasKernGradient(kern, x, covGrad)

% BIASKERNGRADIENT Gradient of bias kernel's parameters.

% IVM

g = sum(sum(covGrad))*kern.variance;