function mlpardKernDisplay(kern)

% MLPARDKERNDISPLAY Display parameters of multi-layer perceptron ARD kernel.

% IVM

fprintf('MLP ARD kernel Variance: %2.4f\n', kern.variance)
fprintf('MLP ARD weight variance: %2.4f\n', kern.weightVariance)
fprintf('MLP ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf('MLP ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
