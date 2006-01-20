function polyardKernDisplay(kern, spacing)

% POLYARDKERNDISPLAY Display parameters of polynomial ARD kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing)
fprintf('Polynomial ARD kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing)
fprintf('Polynomial ARD weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing)
fprintf('Polynomial ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf(spacing)
  fprintf('Polynomial ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
