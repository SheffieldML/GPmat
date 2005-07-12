function linardKernDisplay(kern)

% LINARDKERNDISPLAY Display parameters of linear ARD kernel.

% KERN

fprintf('Linear ARD kernel Variance: %2.4f\n', kern.variance)
for i = 1:kern.inputDimension
  fprintf('Linear ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
