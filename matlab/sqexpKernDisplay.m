function sqexpKernDisplay(kern)

% SQEXPKERNDISPLAY Display parameters of squared exponential kernel.

% KERN

% KERN


fprintf('RBF Variance: %2.4f\n', kern.rbfVariance)
fprintf('RBF inverse width: %2.4f\n', kern.inverseWidth)
fprintf('White noise Variance: %2.4f\n', kern.whiteVariance)
fprintf('Bias Kernel Variance: %2.4f\n', kern.biasVariance)
