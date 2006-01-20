function mlpardKernDisplay(kern, spacing)

% MLPARDKERNDISPLAY Display parameters of multi-layer perceptron ARD kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD kernel Variance: %2.4f\n', kern.variance)
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD weight variance: %2.4f\n', kern.weightVariance)
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  spacing = char(spacing);
  fprintf(spacing);
  fprintf('MLP ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
