function ardKernDisplay(kern)

% ARDKERNDISPLAY Display parameters of ARD kernel.

% KERN


fprintf('RBF Variance: %2.4f\n', kern.rbfVariance)
fprintf('RBF inverse width: %2.4f\n', kern.inverseWidth)
fprintf('Linear Variance: %2.4f\n', kern.linearVariance)
fprintf('Noise Variance: %2.4f\n', kern.whiteVariance)
fprintf('Bias Variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf('Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
