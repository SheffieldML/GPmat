function g = linKernGradient(kern, x, covGrad)

% LINKERNGRADIENT Gradient of lin kernel's parameters.

% IVM

linPart = linKernCompute(kern, x);
g(1) = sum(sum(covGrad.*linPart))/kern.variance;
