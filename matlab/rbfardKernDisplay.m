function rbfardKernDisplay(kern)

% RBFARDKERNDISPLAY Display parameters of radial basis function ARD kernel.

% KERN

% KERN


fprintf('RBF ARD Variance: %2.4f\n', kern.variance)
fprintf('RBF ARD inverse width: %2.4f\n', kern.inverseWidth)
for i = 1:kern.inputDimension
  fprintf('RBF ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
