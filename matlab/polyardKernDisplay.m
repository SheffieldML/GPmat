function polyardKernDisplay(kern)

% POLYARDKERNDISPLAY Display parameters of polynomial ARD kernel.

% KERN

fprintf('Polynomial ARD kernel Variance: %2.4f\n', kern.variance)
fprintf('Polynomial ARD weight variance: %2.4f\n', kern.weightVariance)
fprintf('Polynomial ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf('Polynomial ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
