function linardKernDisplay(kern, spacing)

% LINARDKERNDISPLAY Display parameters of linear ARD kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing)
fprintf('Linear ARD kernel Variance: %2.4f\n', kern.variance)
for i = 1:kern.inputDimension
  fprintf(spacing)
  fprintf('Linear ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
