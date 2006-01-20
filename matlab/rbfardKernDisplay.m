function rbfardKernDisplay(kern, spacing)

% RBFARDKERNDISPLAY Display parameters of radial basis function ARD kernel.

% KERN


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF ARD Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('RBF ARD inverse width: %2.4f\n', kern.inverseWidth)
for i = 1:kern.inputDimension
  fprintf(spacing);
  fprintf('RBF ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
