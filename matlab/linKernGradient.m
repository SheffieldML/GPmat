function g = linKernGradient(kern, x, covGrad)

% LINKERNGRADIENT Gradient of lin kernel's parameters.

% IVM

linPart = linKernCompute(kern, x);
if kern.linearBound
  g(1) = sum(sum(covGrad.*linPart))/kern.variance*gradFactLinearBound(kern.variance);
else
  g(1) = sum(sum(covGrad.*linPart));
end