function g = linKernGradient(kern, x, covGrad)

% LINKERNGRADIENT Gradient of lin kernel's parameters.

% IVM

linPart = linKernCompute(x, kern);
g(1) = sum(sum(covGrad.*linPart));