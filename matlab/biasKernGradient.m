function g = biasKernGradient(kern, x, covGrad)

% BIASKERNGRADIENT Gradient of bias kernel's parameters.

% IVM
if kern.linearBound
  g = sum(sum(covGrad))*gradFactLinearBound(kern.variance);
else
  g = sum(sum(covGrad))*kern.variance;
end