function g = whiteKernGradient(kern, x, covGrad)

% WHITEKERNGRADIENT Gradient of white noise kernel's parameters.

% IVM

g = trace(covGrad);
