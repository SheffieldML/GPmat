function g = whiteKernGradient(kern, x, covGrad)

% WHITEKERNGRADIENT Gradient of white noise kernel's parameters.

% IVM
if kern.linearBound
  g = trace(covGrad)*gradFactLinearBound(kern.variance);
else
  g = trace(covGrad)*kern.variance;
end
